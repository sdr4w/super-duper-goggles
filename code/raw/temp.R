library(raster);
library(data.table);
library(spatstat);
library(imager);
library(akima);
library(compiler);
library(oro.dicom);
library(oro.nifti);

############################
##        Constant        ##
############################

INPUT_FOLDER <- "data/train/sample_images";
TRUTH_LABELS <- fread("./data/train/stage1_labels.csv");

############################
## User-Defined Functions ##
############################


#' @name        fn_Age
#' @author      Steven Rankine, Kendall-Greene Associates
#' @date        2017-02-11 
#' @description  
#' @references   
#' 
fn_Age <- cmpfun(function(birthDate, output="y"){
  t2 <- as.integer(strsplit(format(Sys.time(), "%Y-%m-%d"), "-")[[1]]);
  t1 <- as.integer(strsplit(str2date(birthDate, format.in="%Y%m%d", format.out="%Y-%m-%d"), "-")[[1]]);
  t0 <- abs(t2 - t1);
  t0 <- t0 * switch(output,
    d = c(365.000000000, 30.416670000,  1.000000000),
    m = c( 12.000000000,  1.000000000,  0.002739726),
    y = c(  1.000000000,  0.083333333,  0.002739726)
  );
  return (sum(t0));
});


#' @name        fn_Malignancy
#' @author      Steven Rankine, Kendall-Greene Associates
#' @date        2017-02-11 
#' @description  
#' @references   
#' 
fn_Malignancy <- cmpfun(function(nodeSize){
  x_min <- min(nodeSize);
  x_max <- max(nodeSize);
  mayo  <- 0.1;
  elcap <- 0.0;
  if(x_min > 2  & x_max < 5) {
    mayo  <- 0.1;
    elcap <- 1.0;
  } else if(x_min > 4  & x_max < 7) {
    mayo  <- 0.7;
    elcap <- 1.0;
  } else if(x_min > 6  & x_max < 10) {
    mayo  <- 0.7;
    elcap <- 24.0;
  } else if(x_min > 8  & x_max < 20) {
    mayo  <- 18.7;
    elcap <- 24.0;
  } else if(x_min > 21 & x_max < 30) {
    mayo  <- 33.3;
    elcap <- 24.0;
  } else if(x_min > 21 & x_max < 45) {
    mayo  <- 33.3;
    elcap <- 80.0;
  } else if(x_min > 45) {
    mayo  <- 33.3;
    elcap <- 80.0;
  }
  return (list(elcap=elcap,mayo=mayo));
});


#' @name        fn_Resample
#' @author      Steven Rankine, Kendall-Greene Associates
#' @date        2017-02-03 
#' @description  
#' @references   
#' 
fn_Resample <- cmpfun(function(im, old_spc, new_spc) {
  
  zo <- im;
  if( (old_spc[1] != new_spc[1]) | (old_spc[2] != new_spc[2]) ){
    n  <- dim(im);
    x  <- (1:n[1]) * old_spc[1];
    y  <- (1:n[2]) * old_spc[2];
    z  <- im;
    xo <- (1:n[1]) * new_spc[1];
    yo <- (1:n[2]) * new_spc[2];
    li <- bicubic.grid(x, y, z, dx=new_spc[1], dy=new_spc[2]);
    li <- bicubic.grid(li$x, li$y, li$z, nx=n[1], ny=n[2]);
    zo <- li$z;
  }
  return (zo);
  
});


#' @name        fn_CT2HU
#' @author      Steven Rankine, Kendall-Greene Associates
#' @date        2017-01-31 
#' @description  
#' @references   
#' 
fn_CT2HU <- cmpfun(function(im, slope, intcp) {
  
  slope <- as.double(slope);
  intcp <- as.integer(intcp);
  
  ## Set outside-of-scan pixels to -1000 (air)
  im[im <= -2000] <- -1000;
  im[im >= +2000] <- +1000;
  im <- (im * slope) + intcp;

  return (im);
  
});


#' @name        fn_DetectThreshold
#' @author      Steven Rankine, Kendall-Greene Associates
#' @date        2017-01-31 
#' @description  
#' @references   
#' 
fn_DetectThreshold <- cmpfun(function(im){ 

  ## Filter out non-lung tissue  
  im[im <  -700] <- -1000; # Air
  im[im >= +400] <- +1000; # Bone, Blood, Muscle, etc

  ## Apply a Gaussian blur to reduce noise and avoid false circle detection:
  im <- blur(as.im(im))$v;
  
  ## Statistics of image center
  lb <- as.integer(dim(im)*0.25);
  ub <- as.integer(dim(im)*0.75);
  im_ctr <- as.double(im[lb[1]:ub[1], lb[2]:ub[2]]);
  rm(lb, ub);

  ## Remove underflow bins
  im[im == max(im_ctr) | im == min(im_ctr)] <- mean(im_ctr);
  # im <- (im-mean(im_ctr))/sd(im_ctr);
  
  ## Find clusters
  set.seed(1989);
  im_clust <- kmeans(as.double(im[im != -1000 & im != 1000 & im != 0]), 2);
  im_thres <- mean(im_clust$centers);
  im[im <  im_thres] <- 0;
  im[im >= im_thres] <- 1;
  rm(im_clust, im_thres, im_ctr); 
  
  return (im);
  
});


#' @name        fn_ProcessDicom
#' @author      Steven Rankine, Kendall-Greene Associates
#' @date        2017-01-27 
#' @description Extract header info from DICOM data structure 
#'              1-Convert to numeric where possible
#'              2-Add missing data structures
#' @references   
#' 
fn_ProcessDicom <- cmpfun(function(dcm){
  
  ## 
  hdr    <- dicomTable(dcm$hdr);
  img_ct <- dcm$img;
  img_hu <- list(NULL);   
  img_th <- list(NULL);
  info   <- list(NULL);
  imax   <- length(img_ct);

  ## Convert DICOM fields to numeric where possible
  hdr$`0020-0013-InstanceNumber`      <- as.numeric(hdr$`0020-0013-InstanceNumber`);
  hdr$`0020-1041-SliceLocation`       <- as.numeric(hdr$`0020-1041-SliceLocation`);
  hdr$`0028-0010-Rows`                <- as.integer(hdr$`0028-0010-Rows`);
  hdr$`0028-0011-Columns`             <- as.integer(hdr$`0028-0011-Columns`);
  hdr$`0028-1052-RescaleIntercept`    <- as.double( hdr$`0028-1052-RescaleIntercept`);
  hdr$`0028-1053-RescaleSlope`        <- as.double( hdr$`0028-1053-RescaleSlope`);

  ## Sort into slice order
  if(is.null(hdr$`0020-0013-InstanceNumber`)){
    tmpZ <- hdr$`0020-0032-ImagePositionPatient`;
    tmpZ <- matrix(as.numeric(unlist(strsplit(tmpZ," "))), ncol=3, byrow=T)[,3];
  } else {
    tmpZ <- hdr$`0020-0013-InstanceNumber`;
  }
  hdr <- hdr[order(tmpZ),];      # Sort metadata
  img_ct <- img_ct[order(tmpZ)]; # sort original CT images

  ## Intialize Feature Variables
  info$Lung.Normal   <- rep(0, imax);
  info$Lung.Benign   <- rep(0, imax);
  info$Lung.Malign   <- rep(0, imax);
  info$Patient.Age   <- rep(0, imax);
  info$Slice.Space.x <- rep(0, imax);
  info$Slice.Space.y <- rep(0, imax);
  info$Slice.Space.z <- rep(abs(tmpZ[order(tmpZ)][1]-tmpZ[order(tmpZ)][2]), imax);
  rm(tmpZ);
  
  ## Process images & extract features
  spc2 <- c(0.597656, 0.597656);
  for(i in 1:imax) {
    slope <- as.double( hdr$`0028-1053-RescaleSlope`[i]);
    intcp <- as.integer(hdr$`0028-1052-RescaleIntercept`[i]);
    spc1  <- as.double(unlist(strsplit(hdr$`0028-0030-PixelSpacing`[[1]]," ")));
    ## Save transformations
    # img_ct[[i]] <- fn_Resample(img_ct[[i]], spc1, spc2);       # 
    img_hu[[i]] <- fn_CT2HU(img_ct[[i]], slope, intcp);          # Convert CT scans into Hounsfield units (HU)
    img_th[[i]] <- fn_DetectThreshold(img_hu[[i]]);              # Convert HU scan into B&W
    
    ## 
    circles <- fn_HoughCircle(img_th[[i]],spc1);

    ## HU Pixel Counting Features
    tmp <- img_hu[[i]];
    info$Lung.Normal[i] <- length(tmp[tmp < -300 & tmp > -700]);# Possible Normal Lung Tissue
    info$Lung.Benign[i] <- length(tmp[tmp <    0 & tmp > -20]); # Possible Benign Lung Nodule Tissue
    info$Lung.Malign[i] <- length(tmp[tmp < +200 & tmp >   0]); # Possible Malignant Lung Nodule Tissue
    ## Patient Features
    info$Patient.Age[i]  <- fn_Age(hdr$`0010-0030-PatientsBirthDate`[i]);
    info$Patient.Risk[i] <- FALSE;
    ## 
    tmp <- fn_Malignancy(circles[,3]);
    info$Diagnosis.Mayo[i]  <- tmp$mayo;
    info$Diagnosis.ELCAP[i] <- tmp$elcap;
    ## Slice Feature
    info$Slice.Space.x[i] <- spc1[1];
    info$Slice.Space.y[i] <- spc1[2];
    rm(tmp, circles, spc1, slope, intcp);
  }
  rm(hdr, imax, spc2); gc();
  return(list(info=info,img_ct=img_ct,img_hu=img_hu,img_th=img_th));
  
});


#' @name        fn_IsCenterOfCircle
#' @author      Steven Rankine, Kendall-Greene Associates
#' @date        2017-02-13 
#' @description  
#' @references   
#' 
fn_IsCenterOfCircle <- cmpfun(function(a, b, rmin, rmax, im, cvalue=1) {
  
  im <- as.matrix(im);
  out <- c(a,b,0);
  if(im[a,b] == cvalue){
    tht <- seq(from=0, to=2*pi, by=pi/16);
    # nextR <- TRUE;
    for(R in rmin:rmax){

      x <- as.integer(a+R*cos(tht));
      y <- as.integer(b+R*sin(tht));
      v <- extract(im,cbind(x,y));
      if( (length(v==cvalue)/length(v)) < 0.8){ break; }
      out[3] <- R; 
      
                  
      # for(theta in seq(from=0, to=2*pi, by=pi/8)){
      #   z <- c(as.integer(a+R*cos(theta)), as.integer(b+R*sin(theta)));
      #   z[z<c(1,1)]  <- 1;
      #   z[z>dim(im)] <- min(dim(im));
      #   nextR <- nextR & ifelse(im[z[1],z[2]] == im[a,b], T, F);
      # }
      # if(!nextR){ break; }
      # out[3] <- R; 
    }
  }
  return (out);
  
});

fn_HoughCircle <- cmpfun(function(im, spc){
  circles <- data.table(a=NULL,b=NULL,R=NULL);
  lb  <- as.integer(dim(im)*0.3);
  ub  <- as.integer(dim(im)*0.7);
  rmx <- as.integer(45.0/min(spc));
  rmn <- as.integer(4.0/min(spc));
  for(j in lb[1]:ub[1]){
    for(k in lb[2]:ub[2]){
      tmp <- fn_IsCenterOfCircle(j, k, rmn, rmx, im);
      if(tmp[3] != 0){
        circles <- rbind(circles, tmp);
      }
  }}
  return (circles);
});

k <- 1;
im <- load.image("./code/raw/test.jpg");
im_gray   <- grayscale(im);
im_blur   <- blur_anisotropic(im_gray, ampl=1e10);
im_binary <- round(im_blur/(k*max(im_blur)), 0);
par(mfrow=c(3, 1));
# plot(im, main="Original");
plot(im_gray, main="Greyscale");
plot(im_blur, main="Blur");
plot(im_binary, main="Binary");
# blur_anisotropic(grayscale(im),ampl=1e9,sharp=1) %>% plot(main="Blurred (anisotropic)")
# blur_anisotropic(round(grayscale(im)/max(grayscale(im)),0),ampl=1e9) %>% plot(main="Binary")

par(mfrow=c(3, 1));
round(im_blur/(0.2*max(im_blur)), 0) %>% plot(main="Binary .2");
round(im_blur/(0.233*max(im_blur)), 0) %>% plot(main="Binary .233");
round(im_blur/(0.267*max(im_blur)), 0) %>% plot(main="Binary .267");

a<-300;
b<-100;
R<-75;

theta <- seq(from=0, to=2*pi, by=pi/8);
x  <- as.integer(a+R*cos(theta));
y  <- as.integer(b+R*sin(theta));
im <- as.matrix(im_binary);
tmp  <- cbind(x,y,v=extract(im,cbind(x,y)));
test <- count(tmp[,"v"] == im[a,b]);
summary(test)

############################
##      Main Program      ##
############################


## Create empty data frame to capture scan features
df <- data.frame(
  id      = as.character(NULL), # Unique patient ID
  truth   = as.double(NULL),    # True cancer prognosis
  cancer  = as.double(NULL),    # Predicted cancer prognosis
  ncount  = as.integer(NULL),   # Number of nodules > 3mm
  nwidth  = as.integer(NULL),   # 
  nlength = as.integer(NULL),   # 
  norm    = as.integer(NULL),   # Number of pixels in normal lung range
  benign  = as.integer(NULL),   # Number of pixels in benign nodule range
  malign  = as.integer(NULL),   # Number of pixels in malignant nodule range
  mayo    = as.double(NULL),    # 
  elcap   = as.double(NULL),    # 
  age     = as.double(NULL),    # Patient age as of today
  risk    = as.logical(NULL)    # Is patient a smoker or have family history
);

## Search for features
timestamp();
for(pid in list.files(path=INPUT_FOLDER)){
  dcm   <- readDICOM(paste0(INPUT_FOLDER, "/", pid), verbose=T, recursive=F);
  truth <- TRUTH_LABELS[TRUTH_LABELS$id==pid,'cancer'];
  this  <- fn_ProcessDicom(dcm);
  df <- rbind(df, data.frame(
    id      = as.character(pid),   
    truth   = as.double(ifelse(is.null(nrow(truth)),NA,truth)),  
    cancer  = 0.5,    
    ncount  = 0,     
    nwidth  = "",     
    nlength = "",    
    norm    = as.integer(sum(this$info$Lung.Normal)),     
    benign  = as.integer(sum(this$info$Lung.Benign)),      
    malign  = as.integer(sum(this$info$Lung.Malign)),      
    mayo    = as.double(max(this$info$Diagnosis.Mayo)),   
    elcap   = as.double(max(this$info$Diagnosis.ELCAP)),    
    age     = as.double(unique(this$info$Patient.Age)),  
    risk    = as.logical(unique(this$info$Patient.Risk))    
  ));
  rm(dcm, this, truth);
  gc();
}
timestamp();





# patientID <- list.files(path=INPUT_FOLDER);
# fname <- paste0(INPUT_FOLDER,"/00cba091fa4ad62cc3200a657aeb957e");
# dcm  <- readDICOM(fname, verbose=T, recursive=FALSE);
# this <- fn_ProcessDicom(dcm);
# 
# j <- 1;
# par(mfrow=c(3, 1));
# hist(this$img_ct[[j]]);
# hist(this$img_hu[[j]]);
# hist(this$img_th[[j]]);

par(mfrow=c(3, 1));
image(this$img_ct[[j]], col=gray((0:256)/256)); 
image(this$img_hu[[j]], col=gray((0:256)/256));
image(this$img_th[[j]], col=gray((0:256)/256)); 

