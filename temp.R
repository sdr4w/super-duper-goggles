library(imager);
library(akima);
library(compiler);
library(oro.dicom);
library(oro.nifti);

############################
##        Constant        ##
############################

INPUT_FOLDER <- "data/train/sample_images";

############################
## User-Defined Functions ##
############################

fn_Malignancy <- cmpfun(function(nodeSize){
  x_min <- min(nodeSize);
  x_max <- max(nodeSize);
  if(x_min <= 3){
    mayo  <- 0.1
    elcap <- 0.0;
  } else if(x_min > 2  & x_max < 5) {
    mayo  <- 0.1
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
  }
  return (list(elcap,mayo));
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
  im[im <  -700]              <- -1000; # Air
  # im[im >= -700 & im <  -300] <-  -500; # Normal Lung Tissue
  # im[im >=  -20 & im <     0] <-   +15; # Lung Nodule Tissue
  # im[im >     0 & im <= +200] <-   +15; # Lung Nodule Tissue
  im[im >= +200]              <- +1000; # Bone, Blood, Muscle, etc

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
  today <- as.integer(format(Sys.time(), "%Y%m%d"));
  hdr <- dicomTable(dcm$hdr);
  img_ct <- dcm$img;
  img_hu <- list(NULL);   
  img_th <- list(NULL);
  ## Convert DICOM fields to numeric where possible
  hdr$`0002-0000-GroupLength`         <- as.numeric(hdr$`0002-0000-GroupLength`);
  hdr$`0010-0030-PatientsBirthDate`   <- as.integer(hdr$`0010-0030-PatientsBirthDate`);
  hdr$`0020-0011-SeriesNumber`        <- as.numeric(hdr$`0020-0011-SeriesNumber`);
  hdr$`0020-0012-AcquisitionNumber`   <- as.numeric(hdr$`0020-0012-AcquisitionNumber`);
  hdr$`0020-0013-InstanceNumber`      <- as.numeric(hdr$`0020-0013-InstanceNumber`);
  hdr$`0020-1041-SliceLocation`       <- as.numeric(hdr$`0020-1041-SliceLocation`);
  hdr$`0028-0002-SamplesperPixel`     <- as.numeric(hdr$`0028-0002-SamplesperPixel`);
  hdr$`0028-0010-Rows`                <- as.integer(hdr$`0028-0010-Rows`);
  hdr$`0028-0011-Columns`             <- as.integer(hdr$`0028-0011-Columns`);
  hdr$`0028-0100-BitsAllocated`       <- as.integer(hdr$`0028-0100-BitsAllocated`);
  hdr$`0028-0101-BitsStored`          <- as.integer(hdr$`0028-0101-BitsStored`);
  hdr$`0028-0102-HighBit`             <- as.integer(hdr$`0028-0102-HighBit`);
  hdr$`0028-0103-PixelRepresentation` <- as.integer(hdr$`0028-0103-PixelRepresentation`);
  hdr$`0028-0120-PixelPaddingValue`   <- as.integer(hdr$`0028-0120-PixelPaddingValue`);
  hdr$`0028-1050-WindowCenter`        <- as.double( hdr$`0028-1050-WindowCenter`);
  hdr$`0028-1051-WindowWidth`         <- as.double( hdr$`0028-1051-WindowWidth`);
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

  ## User-defined metadata
  hdr$PatientAge     <- (today - this$hdr$`0010-0030-PatientsBirthDate`)/365;
  hdr$SliceThickness <- abs(tmpZ[order(tmpZ)][1]-tmpZ[order(tmpZ)][2]);
  if(!is.null(hdr$`0020-1041-SliceLocation`)){
    hdr$SliceThickness <- abs(hdr$`0020-1041-SliceLocation`[1] - hdr$`0020-1041-SliceLocation`[2]);
  }
  ## Process images
  spc2 <- c(0.597656, 0.597656);
  for(i in 1:length(img_ct)) {
    slope <- as.double( hdr$`0028-1053-RescaleSlope`[i]);
    intcp <- as.integer(hdr$`0028-1052-RescaleIntercept`[i]);
    spc1  <- as.double(unlist(strsplit(hdr$`0028-0030-PixelSpacing`[[1]]," ")));
    ## Save transformations
    img_ct[[i]] <- fn_Resample(img_ct[[i]], spc1, spc2);# 
    img_hu[[i]] <- fn_CT2HU(img_ct[[i]], slope, intcp); # Convert CT scans into Hounsfield units (HU)
    img_th[[i]] <- fn_DetectThreshold(img_hu[[i]]);     # Convert HU scan into B&W
  }
  
  return(list(hdr=hdr,img_ct=img_ct,img_hu=img_hu,img_th=img_th));
  
});


############################
##      Main Program      ##
############################

df <- data.frame(
  PatientID       = NULL,
  PatientAge      = NULL,
  PatientRisk     = NULL,
  MalignancyELCAP = NULL,
  MalignancyMayo  = NULL
);
timestamp();
for(pid in list.files(path=INPUT_FOLDER)){
  fn   <- paste0(INPUT_FOLDER, "/", pid);
  dcm  <- readDICOM(fn, verbose=T, recursive=FALSE);
  this <- fn_ProcessDicom(dcm);
  df   <- rbind(df,cbind(
    PatientID       = pid,
    PatientAge      = (today - unique(this$hdr$`0010-0030-PatientsBirthDate`))/365,
    PatientRisk     = 0.0,
    MalignancyELCAP = 0.0,
    MalignancyMayo  = 0.0
  ));
}
timestamp();

patientID <- list.files(path=INPUT_FOLDER);
fname <- paste0(INPUT_FOLDER,"/00cba091fa4ad62cc3200a657aeb957e");
dcm  <- readDICOM(fname, verbose=T, recursive=FALSE);
this <- fn_ProcessDicom(dcm);

j <- 1;
par(mfrow=c(3, 1));
hist(this$img_ct[[j]]);
hist(this$img_hu[[j]]);
hist(this$img_th[[j]]);

par(mfrow=c(3, 1));
image(this$img_ct[[j]], col=gray((0:256)/256)); 
image(this$img_hu[[j]], col=gray((0:256)/256));
image(this$img_th[[j]], col=gray((0:256)/256)); 

