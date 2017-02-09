library(zoo);
library(fields);
library(stats);
library(spatstat);
library(oro.dicom);
library(oro.nifti);
library(plotrix);

############################
## User-Defined Functions ##
############################

#' Name...: fn_Diagnosis
#' Author.: Steven Rankine, Systems Engineer
#'          Kendall-Greene Associates
#' Date...: 2016-01-28
#' Purpose: Calculate ejection fraction (EF) from and characterize 
#'          estimate with a diagnosis.
#'          
#' @param  Diastole - Volume during the phase when heart muscle 
#'                    relaxes and allows the chambers to fill 
#'                    with blood.
#' @param  Systole  - Volume during the phase when heart muscle 
#'                    contracts and pumps blood from the chambers 
#'                    into the arteries.
#' @param  Age      - Age of the patient.
#' @param  Sex      - Sex of the patient.
#' @return List that contains the expected normal ejection fraction 
#'         for age/sex (EFNorm), calculated ejection fraction (EF), 
#'         and diagnosis statement. 
#' @export
#' 
fn_Diagnosis <- function(Diastole, Systole, Age, Sex) {
  
  # Calculate Ejection Fraction
  EF <- 100.0 * (Diastole-Systole) / Diastole; 
  
  # Determine descriptive diagnosis, based on range that EF falls into
  iEF <- as.integer(EF); # Round EF to the nearesr integer 
  if(iEF > 75)            { Diagnosis <- "Hyperdynamic"; } 
  else if(iEF %in% 65:74) { Diagnosis <- "Normal to Hyperdynamic"; } 
  else if(iEF %in% 55:64) { Diagnosis <- "Normal EF"; } 
  else if(iEF %in% 45:54) { Diagnosis <- "Mildly Abnormal"; } 
  else if(iEF %in% 35:44) { Diagnosis <- "Moderately Abnormal"; } 
  else if(iEF %in%  0:34) { Diagnosis <- "Seriously Abnormal"; }
  else                    { Diagnosis <- "Dead On Arrival"; }
  
  # Determine the average ejection fraction based on sex and age.  
  # Based on:"Age and gender specific normal values of left ventricular  
  # mass, volume and function for gradient echo magnetic resonance   
  # imaging: a cross sectional study", By Peter A Cain, Ragnhild Ahl,   
  # Erik Hedstrom, Martin Ugander, Ase Allansdotter-Johnsson, Peter   
  # Friberg, Hakan Arheden; BMC Med Imaging. 2009; 9: 2. Published o  
  # nline 2009 January 21. doi: 10.1186/1471-2342-9-2, PMCID: PMC2657902
  
  EFNorm <- 67.0;
  Age <- as.integer(substr(Age,1,nchar(Age)-1));
  if(Sex == "M") {
    if(Age %in% 11:20)      { EFNorm <- 67.0; } 
    else if(Age %in% 11:30) { EFNorm <- 66.0; }
    else if(Age %in% 31:50) { EFNorm <- 65.0; }
    else if(Age %in% 51:70) { EFNorm <- 64.0; }
    else if(Age %in% 71:80) { EFNorm <- 63.0; }
    else if(Age > 80)       { EFNorm <- 62.0; }
  } else if(Sex == "F") {
    if(Age %in% 11:20)      { EFNorm <- 71.0; }
    else if(Age %in% 21:40) { EFNorm <- 70.0; }
    else if(Age %in% 41:70) { EFNorm <- 69.0; }
    else if(Age %in% 71:80) { EFNorm <- 68.0; }
    else if(Age >80)        { EFNorm <- 67.0; } 
  }
  
  return (list(NormalEF=EFNorm, EF=EF, Diagnosis=Diagnosis, Diastole=Diastole, Systole=Systole));
  
}


#' Name...: fn_FrustrumVolume
#' Author.: Steven Rankine, Systems Engineer
#'          Kendall-Greene Associates
#' Date...: 2016-01-27
#' Purpose: This function calculate the Volume in milli-liters
#'           of circular cone frustrum for each segment. 
#' @param df        - The images (data.frame)
#' @param thickness - Height of segment in mm (integer)
#' @return 
#' @export
#' 
fn_FrustrumVolume <- function( df, thickness ) {
  
  result <- array(dim=c(dim(df)[1],dim(df)[2]-1));
  for(i in 1:dim(df)[1]) {
    h <- thickness[i];
    for(j in 1:dim(df)[2]) {
      r2 <- df[i,j];
      if(j != 0){
        r1 <- df[i,j-1];
        result[i,j-1] <- 0.001*pi*h*(r1**2 + r1*r2 + r2*2)/3; # Volume (ml) of circular cone frustrum
      }
    }  
  }  
  rownames(result) <- rownames(result, do.NULL=FALSE, prefix="obs.");
  colnames(result) <- colnames(result, do.NULL=FALSE, prefix="seg.");
  return (result);
  
}


#' Calculate average of a time-series collection of images
#' Author.: Steven Rankine, Systems Engineer
#'          Kendall-Greene Associates
#' Date...: 2016-01-27
#' Purpose: This function takes a time-series collection of images 
#'          and calculates the DC component (mean over time) of 
#'          the pixel intensities.
#' @param imgs    - The data we are querying (data.frame)
#' @param Rows    - Vector of row lengths of each image (integer)
#' @param Columns - Vector of clumns lengths of each image  (integer)
#' @return A data structure of matrices (Average, Variance, and blurred variance)
#' @export
#' 
fn_DCBias <- function(imgs, Rows, Columns) {
  
  # Loop through time-series collection of images to create a single 
  # 3-diminsional matrix of pixel intensity values
  tmp <- array(dim=c(max(Rows),max(Columns),length(imgs)));
  for(k in 1:length(imgs)) {
    for(i in 1:Rows[k]) {
      for(j in 1:Columns[k]) {
        tmp[i,j,k] <- (imgs[[k]])[i,j];
      }
    }
  }
  
  # Return value
  result <- NULL;                                            # Initialize return value
  result$avg <- apply(tmp, c(1,2), mean, na.rm=TRUE);        # Matrix of Average over time (DC Comp/Bias/Offset)
  result$var <- apply(tmp, c(1,2), var,  na.rm=TRUE);        # Matrix of Variance over time
  result$roi <- (blur(as.im(result$var), normalise=TRUE))$v; # Matrix of Gaussian blur of variance
  
  # Determine center
  result$tst <- result$roi/max(result$roi);
  result$tst[result$tst <= 0.250] <- 0.0;
  result$tst[result$tst >  0.250 && result$tst < 0.375] <- 0.25;
  result$tst[result$tst >= 0.375 && result$tst < 0.625] <- 0.50;
  result$tst[result$tst >= 0.625 && result$tst < 0.875] <- 0.75;
  result$tst[result$tst >= 0.875] <- 1.0;
  
  result$ellipse <- list(
    x=as.integer(0.5*dim(result$tst)[1]), 
    y=as.integer(0.5*dim(result$tst)[2]), 
    a=1, b=1
  );
  
  for(i in 1:dim(result$tst)[1]) {
    test <- result$tst[i,];
    filt <- test > 0.25;
    if(!is.na(test[filt]) && length(test[filt]) > result$ellipse$a) {
      result$ellipse$a <- as.integer(0.5*length(test[filt]));
      result$ellipse$y <- which(filt)[1] + result$ellipse$a;
    }
  }
  
  for(i in 1:dim(result$tst)[2]) {
    test <- result$tst[,i];
    filt <- test > 0.25;
    if(!is.na(test[filt]) && length(test[filt]) > result$ellipse$b) {
      result$ellipse$b <- as.integer(0.5*length(test[filt]));
      result$ellipse$x <- which(filt)[1] + result$ellipse$b;
    }
  }
  
  return (result);
  
}


#' Name...: fn_AlignDataFrames
#' Author.: Steven Rankine, Systems Engineer
#'          Kendall-Greene Associates
#' Date...: 2016-01-27
#' Purpose: 
#'           
#'          
#' @param df       - 
#' @param seg      - 
#' @param instance - 
#' @return 
#' @export
#' 
fn_AlignDataFrames <- function( df, seg, instance ) {
  
  instance <- as.character(instance);
  for(i in 1:length(instance)) {
    if(as.integer(instance[i]) < 10) {
      instance[i] <- paste0("0",instance[i]);
    }
  }
  
  rownames(seg) <- instance;
  if(is.null(df)) {
    result <- cbind(df, seg);
  } else {
    result <- cbind(df,rep(NA,each=length(df)));
    j <- dim(result)[2];
    for(i in rownames(seg)) {
      if(!(i %in% rownames(result))) {
        result <- rbind(result,rep(NA,each=j));
        rownames(result)[rownames(result)==""] <- i;
      }
      result[as.integer(i),j] <- seg[as.integer(i)];
    }
  }
  result[order(row.names(result)),];
  return (result);
  
}


#' 
#' Author.: Steven Rankine, Systems Engineer
#'          Kendall-Greene Associates
#' Date...: 2016-01-27
#' Purpose: 
#' @param imgs  - 
#' @param nrows - 
#' @param ncols - 
#' @return 
#' @export
#' 
fn_ZeroFill <- function(imgs, nrows, ncols) {
  
  for(i in 1:length(imgs)) {
    img <- imgs[[i]];
    if((nrows - dim(img)[1]) > 0) {
      delta <- max(1,as.integer((nrows-dim(img)[1])/2.0));
      img <- rbind(matrix(data=0, nrow=delta, ncol=dim(img)[2]), img);
      img <- rbind(img, matrix(data=0, nrow=delta, ncol=dim(img)[2]));
    }
    if((ncols - dim(img)[2]) > 0) {
      delta <- max(1,as.integer((ncols - dim(img)[2])/2.0));
      img <- cbind(matrix(data=0, nrow=dim(img)[1], ncol=delta), img);
      img <- cbind(img, matrix(data=0, nrow=dim(img)[1], ncol=delta));
    }
    imgs[[i]] <- img;
  }
  return (imgs);
  
}


#' Name...: fn_DetectRadius
#' Author.: Steven Rankine, Systems Engineer
#'          Kendall-Greene Associates
#' Date...: 2016-01-27
#' Purpose: This function takes a 2-D image and attempts to detect the
#'          radius of the ventricular structure  
#' Detect the radius of ventricular heart structure from a 2-D image 
#' of a heart cross section
#' @param tmp     - The image (data.frame)
#' @param nrows   - Number of row in image (integer)
#' @param ncols   - Number of columns in image  (integer)
#' @param spacing - Vector of ??? (integer)
#' @return A vector of possible radius values
#' @export
#' 
#fn_DetectRadius <- function(tmp, nrows, ncols, spacing) {
fn_DetectRadius <- function(tmp, nrows, ncols, rowc, colc, spacing) {
    
  result <- rep(0.0, 8);       # Intial vector of return values (Right, Left, Up, and down from center)

  #if(nrows == rowc)
  Q1 <- tmp[rowc:nrows, colc:ncols];
  Q2 <- tmp[rowc:nrows, 1:colc];
  Q3 <- tmp[1:rowc, 1:colc];
  Q4 <- tmp[1:rowc, colc:ncols];
  
  for(i in 1:length(result)) {
    h <- switch(i, 
      spacing[2], sqrt(spacing[1]^2+spacing[2]^2), 
      spacing[2], sqrt(spacing[1]^2+spacing[2]^2), 
      spacing[1], sqrt(spacing[1]^2+spacing[2]^2), 
      spacing[1], sqrt(spacing[1]^2+spacing[2]^2)
    );
    x <- switch(i, 
      tmp[rowc, colc:ncols],  diag(Q1[,seq(ncol(Q1),1)]), 
      rev(tmp[rowc, 1:colc]), diag(Q2), 
      tmp[rowc:nrows, colc],  diag(Q3[,seq(ncol(Q3),1)]), 
      rev(tmp[1:rowc, colc]), diag(Q4) 
    );
    
    x <- x[x!=0.0];
    x[x<(mean(x)-0.25*sd(x))] <- 0.0;
    for(intesity in x) {
      if(intesity > 0.0) {
        result[i] <- result[i] + h;
      } else if(result[i] > 0.0) {
        break;
      }
    }
  }
  result <- result[result!=0];
  return (result);
  
}


fn_Boundary <- function(img) {
  x <- colSums(img);
  y <- rowSums(img);
  return(list(
    xbeg = which(x!=0)[1],
    xend = which(x!=0)[length(which(x!=0))],
    yend = which(y!=0)[length(which(y!=0))],
    ybeg = which(y!=0)[1],
    ctr  = c(
      as.integer(which(x!=0)[1]+(which(x!=0)[length(which(x!=0))]-which(x!=0)[1])*0.5),
      as.integer(which(y!=0)[1]+(which(y!=0)[length(which(y!=0))]-which(y!=0)[1])*0.5)
    )
  ));  
}


fn_SearchPattern <- function(img, xbeg, xend, ybeg, yend, spacing) {
  
  firstPass <- TRUE;

  xbeg <- as.integer(xbeg); 
  xend <- as.integer(xend);
  ybeg <- as.integer(ybeg); 
  yend <- as.integer(yend);
  
  xrng <- min(dim(img)[2],xend-xbeg);
  yrng <- min(dim(img)[1],yend-ybeg);
  
  xsamp <- sample(xbeg:xend, size=xrng);
  ysamp <- sample(ybeg:yend, size=yrng);

  for(x in xsamp[xsamp!=xrng]) {
    for(y in ysamp[ysamp!=yrng]) {
      r <- fn_DetectRadius(img, xrng, yrng, x, y, spacing);
      if(length(r)>1) {
        if(firstPass) {
          firstPass <- FALSE;
          r_var <- var(r);
          result <- c(x, y, mean(r));
        } else {
          if(var(r) < r_var) {
            result <- c(x, y, mean(r));
            r_var  <- var(r);
          }
        }
      }
    }
  }
  return (result);
  
}




#' 
#' Author.: Steven Rankine, Systems Engineer
#'          Kendall-Greene Associates
#' Date...: 2016-01-27
#' Purpose: 
#' @param imgs  - 
#' @param nrows - 
#' @param ncols - 
#' @return 
#' @export
#' 
fn_SearchPattern2 <- function(img) {
  
  r1 <- seq(1, dim(img)[1], dim(img)[1]/16);
  r2 <- seq(1, dim(img)[2], dim(img)[2]/16);

  result <- c(
    as.integer(0.5*dim(img)[1]), 
    as.integer(0.5*dim(img)[2]), 
    mean(img[r1[8:9],r2[8:9]])
  );
  for(i in 1:16){
    tmp <- switch(i, 
      list(img=img[r1[5:6],  r2[5:6]],   ctr=c(mean(r1[5:6]),  mean(r2[5:6]))  ), 
      list(img=img[r1[5:6],  r2[7:8]],   ctr=c(mean(r1[5:6]),  mean(r2[7:8]))  ), 
      list(img=img[r1[5:6],  r2[9:10]],  ctr=c(mean(r1[5:6]),  mean(r2[9:10])) ), 
      list(img=img[r1[5:6],  r2[11:12]], ctr=c(mean(r1[5:6]),  mean(r2[11:12]))), 
      list(img=img[r1[7:8],  r2[5:6]],   ctr=c(mean(r1[7:8]),  mean(r2[5:6]))  ),
      list(img=img[r1[7:8],  r2[7:8]],   ctr=c(mean(r1[7:8]),  mean(r2[7:8]))  ),
      list(img=img[r1[7:8],  r2[9:10]],  ctr=c(mean(r1[7:8]),  mean(r2[9:10])) ),
      list(img=img[r1[7:8],  r2[11:12]], ctr=c(mean(r1[7:8]),  mean(r2[11:12]))), 
      list(img=img[r1[9:10], r2[5:6]],   ctr=c(mean(r1[9:10]), mean(r2[5:6]))  ),
      list(img=img[r1[9:10], r2[7:8]],   ctr=c(mean(r1[9:10]), mean(r2[7:8]))  ),
      list(img=img[r1[9:10], r2[9:10]],  ctr=c(mean(r1[9:10]), mean(r2[9:10])) ), 
      list(img=img[r1[9:10], r2[11:12]], ctr=c(mean(r1[9:10]), mean(r2[11:12]))), 
      list(img=img[r1[11:12],r2[5:6]],   ctr=c(mean(r1[11:12]),mean(r2[5:6]))  ),
      list(img=img[r1[11:12],r2[7:8]],   ctr=c(mean(r1[11:12]),mean(r2[7:8]))  ),
      list(img=img[r1[11:12],r2[9:10]],  ctr=c(mean(r1[11:12]),mean(r2[9:10])) ),
      list(img=img[r1[11:12],r2[11:12]], ctr=c(mean(r1[11:12]),mean(r2[11:12]))) 
    );
    
    if(result[3] < mean(tmp$img)) {
      result <- c(
        as.integer(tmp$ctr[1]), 
        as.integer(tmp$ctr[2]), 
        mean(tmp$img)
      );
    }

  }
  
  return (result);
  
}


##################
## Main Program ##
##################

actual <- read.csv("../data/train.csv");  # read csv file 
final  <- NULL;
#for( PatientID in c(200,346,401) ) { 
for( PatientID in sample(1:500,50) ) { 
  
  baseDir <- paste0("../data/train/",PatientID,"/study");
  short_axial_stack_size <- 0;
  number_of_images       <- 0;
  
  df0 <- NULL;
  df0_colnames <- NULL;
  thick <- NULL;

  par(mfrow=c(2, 1));
  print(baseDir);
  for( i in 1:99 ) { 
    sax <- paste0("sax_",i);
    if( sax %in% list.files(path=baseDir, pattern="sax_*") ) {
      
      # Copy files to temp directory
      dir.create(paste0(baseDir,"/",sax,"/temp"));
      for( fn in list.files(path=paste0(baseDir,"/",sax)) ) {
        # Ignore zero length files
        if(file.info(paste0(baseDir,"/",sax,"/",fn))$size>0) {
          file.copy(from=paste0(baseDir,"/",sax,"/",fn), to=paste0(baseDir,"/",sax,"/temp/",fn));
          number_of_images = number_of_images + 1;
        }
      }
      
      # Read all the DICOM files in temp directory
      dcm <- readDICOM(paste0(baseDir,"/",sax,"/temp"), verbose=FALSE, recursive=FALSE);
      dcm.info <- dicomTable(dcm$hdr);
      unlink((paste0(baseDir,"/",sax,"/temp")), recursive=TRUE);
      
      # PATIENT ATTRIBUTES (C.7-1)
      Patient.Name <- unique(dcm.info["0010-0010-PatientsName"][[1]]);               # Patient's full name.
      Patient.Sex  <- unique(dcm.info["0010-0040-PatientsSex"][[1]]);                # Sex of the named patient. Enumerated Values: M = male F = female O = other
      Patient.Age  <- unique(dcm.info["0010-1010-PatientsAge"][[1]]);                # Age of the Patient.
      
      # GENERAL IMAGE MODULE ATTRIBUTES (C.7-9)
      InstanceNumber   <- as.integer(  dcm.info["0020-0013-InstanceNumber"][[1]]);         # A number that identifies this image. 

      # IMAGE PLANE ATTRIBUTES (C.7-10)
      PixelSpacing     <- as.character(dcm.info["0028-0030-PixelSpacing"][[1]]);           # Physical distance in the patient between the center of each pixel (mm)
      SliceThickness   <- as.integer(  dcm.info["0018-0050-SliceThickness"][[1]]);         # Nominal slice thickness, in mm.
      SeriesDescription<- as.character(dcm.info["0008-103E-SeriesDescription"][[1]]);      # 
      TriggerTime      <- as.numeric(  dcm.info["0018-1060-TriggerTime"][[1]]);            # 

      #ImagePosition    <- as.integer(  dcm.info["0020-0032-ImagePositionPatient"][[1]]);   # The x, y, and z coordinates of the upper left hand corner (center of the first voxel transmitted) of the image, in mm.
      #ImageOrientation <- as.integer(  dcm.info["0020-0037-ImageOrientationPatient"][[1]]);# The direction cosines of the first row and the first column with respect to the patient.
      #SliceLocation    <- as.integer(  dcm.info["0020-1041-SliceLocation"][[1]]);          # Relative position of exposure expressed in mm.
      #StudyTime        <- as.integer(  dcm.info["0008-0030-StudyTime"][[1]]);              # Time the Study started.
      #SeriesTime       <- as.integer(  dcm.info["0008-0031-SeriesTime"][[1]]);             # Time the Series started.
      #RepetitionTime   <- as.integer(  dcm.info["0018-0080-RepetitionTime"][[1]]);         # The period of time in msec between the beginning of a pulse sequence and the beginning of the succeeding pulse sequence
      #NumberOfImages   <- as.integer(  dcm.info["0018-1090-CardiacNumberOfImages"][[1]]);  # 
      
      # 
      Rows             <- as.integer(  dcm.info["0028-0010-Rows"][[1]]);                   # Number of rows in the image
      Columns          <- as.integer(  dcm.info["0028-0011-Columns"][[1]]);                # Number of columns in the image
      SmallestPixel    <- as.integer(  dcm.info["0028-0106-SmallestImagePixelValue"][[1]]);# 
      LargestPixel     <- as.integer(  dcm.info["0028-0107-LargestImagePixelValue"][[1]]); # 

      # Use only Sagital images only
      filterSax <- SeriesDescription == "sax";
      if(TRUE %in% filterSax) {

        dcm$img        <- dcm$img[filterSax];
        Rows           <- Rows[filterSax];
        Columns        <- Columns[filterSax];
        
        PixelSpacing   <- PixelSpacing[filterSax];
        TriggerTime    <- TriggerTime[filterSax];
        SliceThickness <- SliceThickness[filterSax];
        LargestPixel   <- LargestPixel[filterSax];
        InstanceNumber <- InstanceNumber[filterSax];

        # Treat images as a periodic function in the time domain
        short_axial_stack_size <- short_axial_stack_size + 1;
        img_dc <- fn_DCBias(dcm$img, Rows, Columns);                         # Calculate DC Component and variance
        ROIFilter1 <- img_dc$roi < (mean(img_dc$roi) + 1.5*sd(img_dc$roi));  # Filter based on blurred variance values

        # Loop through images of cross-section over time
        sliceRad <- matrix(rep(NA,each=length(dcm$img)));
        for(index in 1:length(dcm$img)) {
          
          # Get the next Image & center point
          tmp <- dcm$img[[index]];    
          spacing <- as.numeric(strsplit(PixelSpacing[index]," ")[[1]]);
          
          # Identify Region og Interest (ROI) & Collect Statistics inside the ROI
          tmp[ROIFilter1] <- 0;                         # Filter based on variance of intensity
          tmp.avg <- mean(tmp[tmp!=0]);                 # Average value inside ROI
          tmp.sd  <- sd(tmp[tmp!=0]);                   # Std Deviation value inside ROI
          tmp.max <- max(tmp[tmp!=0]);                  # Maximum value inside ROI
          tmp.min <- min(tmp[tmp!=0]);                  # Minimum value inside ROI

          if(is.nan(tmp.avg)) tmp.avg <- 0;
          if(is.na(tmp.sd))   tmp.sd  <- 0;
          
          thick <- cbind(thick, SliceThickness);
          bound <- fn_Boundary(tmp);
          
          # Estimate inner ventricular radius at this cross-section plane
          optimal <- fn_SearchPattern2(tmp);
          #print(sprintf("%d: %d, %d, %f",index,optimal[1],optimal[2],optimal[3]));
          rBeg <- (optimal[1])-1;
          cBeg <- (optimal[2])-1;
          rEnd <- (optimal[1])+1;
          cEnd <- (optimal[2])+1;
          r3 <- 0.5*mean(c((rEnd-rBeg)*spacing[1],
                           (cEnd-cBeg)*spacing[2]));    # Estimate radius
          slice    <- tmp[rBeg:rEnd,cBeg:cEnd];         # Initial circular area
          sliceAvg <- mean(slice);                      # Initial average intensity for circular area
          #print(sprintf("%d:  %f, %f, %f",index,sliceAvg,tmp.avg,tmp.sd));
          while( sliceAvg > tmp.avg - (0.5*tmp.sd) ) {  # Loop until falls below some threshold
            # Expand circle bounduary
            rBeg <- rBeg - 1;                           
            rEnd <- rEnd + 1;   
            cBeg <- cBeg - 1;
            cEnd <- cEnd + 1;
            
            
            # Check if circle is not too big
            if(rBeg %in% 1:Rows[index] && cBeg %in% 1:Columns[index] && 
               rEnd %in% 1:Rows[index] && cEnd %in% 1:Columns[index]) {
              sliceAvg <- mean(tmp[rBeg:rEnd,cBeg:cEnd]); #   Calculate new average intensity
              r3 <- 0.5*mean(c((rEnd-rBeg)*spacing[1],    #   Estimate inner radius
                               (cEnd-cBeg)*spacing[2]));
            } else {
              break;
            }
          }
          sliceRad[index] <- r3;
          
          # Estimate outer ventricular radius at this cross-section plane
          #rBeg <- as.integer(Rows[index]/4.0);
          #cBeg <- as.integer(Columns[index]/4.0);
          #rEnd <- Rows[index]-rBeg;
          #cEnd <- Columns[index]-cBeg;
          #optimal <- fn_SearchPattern(tmp, bound$xbeg, bound$xend, bound$ybeg, bound$yend, spacing);
          #sliceRad[index] <- (mean(optimal[3])/min(spacing))/min(Columns[index],Rows[index]);
          #xctr  <- (optimal[1])/Columns[index];
          #yctr  <- (optimal[2])/Rows[index];
          
          
          #print(sprintf("%d -- %f %f %f", index, xctr, yctr, sliceRad[index]));
          
          # Plots
          image(tmp, main=sprintf("%d [%d x %d]",InstanceNumber[index],Rows[index],Columns[index]), col=gray((0:64)/64));          
          abline(h=(0:4)/4.0, col="white", lty="dotted");
          abline(v=(0:4)/4.0, col="white", lty="dotted");
          draw.circle(optimal[1]/Rows[index],optimal[2]/Columns[index], (r3/spacing[1])*(1/256), border="yellow");
          #draw.circle(xctr, yctr, sliceRad[index], border="tan2");
          #draw.ellipse(c(img_dc$ellipse$x/Columns[index]), c(img_dc$ellipse$y/Rows[index]), c(img_dc$ellipse$a/Columns[index]), c(img_dc$ellipse$b/Rows[index]), border="red");
          #drape.plot(1:Rows[index], 1:Columns[index], tmp, border=NA, theta=0, phi=45, main="Spectrogram");
        }
        title(paste0(sax," -- Image #",index), outer=TRUE);
        
        #df0 <- cbind(df0,sliceRad);
        df0 <- fn_AlignDataFrames(df0, sliceRad, InstanceNumber);
        colnames(df0)[dim(df0)[2]] <- sax;
        df0_colnames <- c(df0_colnames,sax);        
      }
      print(sprintf("+ Processed: %6s ... [%d of %d Files] ... Radius ~ %.2f mm", sax, length(dcm$img[filterSax]), length(filterSax), sliceRad[index]));
    }
  }

  # Process collection of radii & calculate volume for each section
  df0[df0==0] <- NA;
  #na.approx(df0, 1:dim(df0)[1]); 
  #df0 <- t(apply(df0,1,na.approx));
  
  colnames(df0) <- df0_colnames;
  rownames(df0) <- rownames(df0, do.NULL=FALSE, prefix="obs.");
  df1 <- fn_FrustrumVolume(df0, thick);
  
  # Display final analysis for this patient
  Zcalc <- fn_Diagnosis(max(rowSums(df1), na.rm=FALSE), min(rowSums(df1), na.rm=FALSE),Patient.Age,Patient.Sex);
  Ztrue <- fn_Diagnosis(actual[actual$Id==PatientID,3],actual[actual$Id==PatientID,2],Patient.Age,Patient.Sex);
  print(sprintf("ANALYSIS | Patient #%s", PatientID));  
  print(sprintf("-----------------------------------------------"));  
  print(sprintf("  Short Axis Stack Size .: %d",             short_axial_stack_size));  
  print(sprintf("  Images Processed ......: %d",             number_of_images));  
  print(sprintf("  Patient Name ..........: %s",             Patient.Name));  
  print(sprintf("  Patient Sex ...........: %s",             Patient.Sex));  
  print(sprintf("  Patient Age ...........: %s",             Patient.Age));  
  print(sprintf("  Diastole Volume .......: %.1f (%.1f) ml", Zcalc$Diastole,  Ztrue$Diastole));
  print(sprintf("  Systole Volume ........: %.1f (%.1f) ml", Zcalc$Systole,   Ztrue$Systole));
  print(sprintf("  Ejection Fraction"));
  print(sprintf("    Normal ..............: %.1f %%",        Zcalc$NormalEF));
  print(sprintf("    Calculated ..........: %.1f (%.1f) %%", Zcalc$EF,        Ztrue$EF));
  print(sprintf("    Diagnosis ...........: %s (%s)",        Zcalc$Diagnosis, Ztrue$Diagnosis));
  print(sprintf("-----------------------------------------------"));
  print(" ");
  final <- rbind(final,cbind(
    PatientID, Patient.Name, Patient.Sex, Patient.Age, 
    Ztrue$Diastole, Ztrue$Systole, Ztrue$EF, Ztrue$Diagnosis, 
    Zcalc$Diastole, Zcalc$Systole, Zcalc$EF, Zcalc$Diagnosis));
  
}

#8-27-34-3-19-10
