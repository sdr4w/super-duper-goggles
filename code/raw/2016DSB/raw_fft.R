library(zoo);
library(stats);
library(spatstat);
library(oro.dicom);
library(oro.nifti);

fn_DCComponent <- function(imgs,Rows,Columns) {
  img_dc <- matrix(0, nrow=max(Rows), ncol=max(Columns));
  for(k in 1:length(imgs)) {
    for(i in 1:Rows[k]) {
      for(j in 1:Columns[k]) {
        print(sprintf("i,j,k --> %d,%d,%d",i,j,k));
        img_dc[i,j] <- img_dc[i,j] + (dcm$img[[k]])[i,j];
      }
    }
  }
  img_dc <- as.integer(img_dc/length(imgs));
  return (img_dc);
}


fn_CvtFFT <- function(cs, sample.rate=1, harmonic=NULL) {
  n  <- length(cs);
  cs <- fft(cs) / n; # normalize
  fn_distance_center <- function(c)signif( Mod(c),        4);
  fn_angle           <- function(c)signif( 180*Arg(c)/pi, 3);  
  df <- data.frame(
    cycle    = 0:(n-1),
    freq     = 0:(n-1) * sample.rate / n,
    strength = sapply(cs, fn_distance_center),
    delay    = sapply(cs, fn_angle)
  );
  if(is.null(harmonic)) {
    value <- sum(df$strength * sin(df$freq*df$cycle+df$delay));
  } else {
    value <- df$strength[harmonic] * sin(df$freq[harmonic]*df$cycle[harmonic]+df$delay[harmonic]);
  }
  #value <- fft(value, inverse=TRUE);
  return (value);
}


fn_GetH1 <- function(imgs,i) {
  ff <- mvfft(imgs[i]);
  h1 <- ff[1];
  result <- np.absolute(fft(h1, inverse=TRUE));
  result <- blur(result);
  #image(blur(as.im(dcm$img[[6]])))
  return(result);  
}


fn_A0 <- function(imgs, Rows, Columns) {
  tmp  <- array(dim=c(max(Rows),max(Columns),length(imgs)));
  for(k in 1:length(imgs)) {
    for(i in 1:Rows[k]) {
      for(j in 1:Columns[k]) {
        tmp[i,j,k] <- (imgs[[k]])[i,j];
      }
    }
  }
  out <- NULL;
  out$avg <- apply(tmp, c(1,2), mean,      na.rm=TRUE);
  #  out$med <- apply(tmp, c(1,2), median,    na.rm=TRUE);
  #  out$sd  <- apply(tmp, c(1,2), sd,        na.rm=TRUE);
  out$var <- apply(tmp, c(1,2), var,       na.rm=TRUE);
  out$roi <- (blur(as.im(out$var),3))$v;
  #  out$fft <- apply(tmp, c(1,2), fn_CvtFFT, sample.rate=1);
  return (out);
}

fn_Harmonize <- function(imgs, n) {
  
  nimg <- length(imgs);
  imgs <- fn_ZeroFill(imgs, n, n);
  
  for(i in 1:nimg) {
    imgs[[i]] <- as.im(imgs[[i]]);
  }
  
  #dcm$img <- fn_CropImage(dcm$img, 100);
  #dcm$img <- fn_ShiftIntensity(dcm$img,max(LargestPixel));
  if(nimg == 30) {
    imgs <- harmonise.im(
      imgs[[1]], imgs[[2]], imgs[[3]], imgs[[4]], imgs[[5]], imgs[[6]], imgs[[7]], imgs[[8]], imgs[[9]], imgs[[10]],
      imgs[[11]],imgs[[12]],imgs[[13]],imgs[[14]],imgs[[15]],imgs[[16]],imgs[[17]],imgs[[18]],imgs[[19]],imgs[[20]],
      imgs[[21]],imgs[[22]],imgs[[23]],imgs[[24]],imgs[[25]],imgs[[26]],imgs[[27]],imgs[[28]],imgs[[29]],imgs[[30]]
    );
  } else if(nimg == 29){
    imgs <- harmonise.im(
      imgs[[1]], imgs[[2]], imgs[[3]], imgs[[4]], imgs[[5]], imgs[[6]], imgs[[7]], imgs[[8]], imgs[[9]], imgs[[10]],
      imgs[[11]],imgs[[12]],imgs[[13]],imgs[[14]],imgs[[15]],imgs[[16]],imgs[[17]],imgs[[18]],imgs[[19]],imgs[[20]],
      imgs[[21]],imgs[[22]],imgs[[23]],imgs[[24]],imgs[[25]],imgs[[26]],imgs[[27]],imgs[[28]],imgs[[29]]
    );
  } else if(nimg == 128){
    imgs <- harmonise.im(
      imgs[[1]], imgs[[2]], imgs[[3]], imgs[[4]], imgs[[5]], imgs[[6]], imgs[[7]], imgs[[8]], imgs[[9]], imgs[[10]],
      imgs[[11]],imgs[[12]],imgs[[13]],imgs[[14]],imgs[[15]],imgs[[16]],imgs[[17]],imgs[[18]],imgs[[19]],imgs[[20]],
      imgs[[21]],imgs[[22]],imgs[[23]],imgs[[24]],imgs[[25]],imgs[[26]],imgs[[27]],imgs[[28]]
    );
  } 
  
  for(i in 1:nimg) {
    imgs[[i]] <- imgs[[i]]$v;
  }
  
  return (imgs);
}


fn_ZoomFilter <- function(img,threshold) {
  for(i in 1:dim(img)[1]) {
    rpct <- length((img[i,])[img[i,]==0])/length(img[i,]);
    for(j in 1:dim(img)[2]) {
      cpct <- length((img[,j])[img[,j]==0])/length(img[,j]);
      if(cpct > threshold && rpct > threshold ) img[i,j] <- 0;
    }
  }
  return (img);
}




final <- NULL;
for( PatientID in 1:1 ) { 
  baseDir <- paste0("../data/train/",PatientID,"/study");
  short_axial_stack_size <- 0;
  number_of_images       <- 0;
  df0 <- NULL;
  df0_colnames <- NULL;
  
  df  <- NULL;
  dfi <- 1;
  print(baseDir);
  for( i in 1:99 ) { 
    sax <- paste0("sax_",i);
    if( sax %in% list.files(path=baseDir, pattern="sax_*") ) {
      print(sprintf("+ Processing Segment: %s ....", sax));
      
      # Copy files to temp directory
      dfCol <- 0;
      short_axial_stack_size = short_axial_stack_size + 1;
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
      PatientsName     <- unique(dcm.info["0010-0010-PatientsName"][[1]]);               # Patient's full name.
      PatientsSex      <- unique(dcm.info["0010-0040-PatientsSex"][[1]]);                # Sex of the named patient. Enumerated Values: M = male F = female O = other
      PatientsAge      <- unique(dcm.info["0010-1010-PatientsAge"][[1]]);                # Age of the Patient.
      
      # GENERAL IMAGE MODULE ATTRIBUTES (C.7-9)
      InstanceNumber   <- as.integer(dcm.info["0020-0013-InstanceNumber"][[1]]);         # A number that identifies this image. 
      
      # IMAGE PLANE ATTRIBUTES (C.7-10)
      PixelSpacing     <- dcm.info["0028-0030-PixelSpacing"][[1]];                       # Physical distance in the patient between the center of each pixel (mm)
      ImagePosition    <- as.integer(dcm.info["0020-0032-ImagePositionPatient"][[1]]);   # The x, y, and z coordinates of the upper left hand corner (center of the first voxel transmitted) of the image, in mm.
      ImageOrientation <- as.integer(dcm.info["0020-0037-ImageOrientationPatient"][[1]]);# The direction cosines of the first row and the first column with respect to the patient.
      SliceThickness   <- as.integer(dcm.info["0018-0050-SliceThickness"][[1]]);         # Nominal slice thickness, in mm.
      SliceLocation    <- as.integer(dcm.info["0020-1041-SliceLocation"][[1]]);          # Relative position of exposure expressed in mm.
      
      StudyTime        <- as.integer(dcm.info["0008-0030-StudyTime"][[1]]);              # Time the Study started.
      SeriesTime       <- as.integer(dcm.info["0008-0031-SeriesTime"][[1]]);             # Time the Series started.
      RepetitionTime   <- as.integer(dcm.info["0018-0080-RepetitionTime"][[1]]);         # The period of time in msec between the beginning of a pulse sequence and the beginning of the succeeding pulse sequence
      NumberOfImages   <- as.integer(dcm.info["0018-1090-CardiacNumberOfImages"][[1]]);  # 
      
      # 
      Rows             <- as.integer(dcm.info["0028-0010-Rows"][[1]]);                   # Number of rows in the image
      Columns          <- as.integer(dcm.info["0028-0011-Columns"][[1]]);                # Number of columns in the image
      SmallestPixel    <- as.integer(dcm.info["0028-0106-SmallestImagePixelValue"][[1]]);# 
      LargestPixel     <- as.integer(dcm.info["0028-0107-LargestImagePixelValue"][[1]]); # 
      
      # 
      print(sprintf("  => Calculating statistics "));
      img_dc <- fn_A0(dcm$img, Rows, Columns);
      ROIFilter <- img_dc$roi < (mean(img_dc$roi) + 2*sd(img_dc$roi));
      #threshold <- (mean((blur(as.im(img_dc$var)))$v) + 2*sd((blur(as.im(img_dc$var)))$v));
      #roiFilter <- (blur(as.im(img_dc$var)))$v < threshold;
      
      # Loop through images of cross-section over time
      print(sprintf("  => Estimate radii"));
      sliceRad <- matrix(rep(NA,each=length(dcm$img)));
      for(index in 1:length(dcm$img)) {    
        # Get the next Image & center point
        tmp  <- dcm$img[[index]];    
        rBeg <- (Rows[index]/2);
        cBeg <- (Columns[index]/2);
        rEnd <- (Rows[index]/2);
        cEnd <- (Columns[index]/2);
        spacing <- as.numeric(strsplit(PixelSpacing[index]," ")[[1]]);
        
        # Image transformation
        tmp[ROIFilter] <- 0;       # Filter based on variance
        tmpAvg <- mean(tmp);
        tmpSD  <- sd(tmp);
        tmp[tmp<tmpAvg] <- 0;
        tmp[tmp>tmpAvg && tmp<(tmpAvg+tmpSD)] <- tmpAvg+tmpSD;
        tmp[tmp>(tmpAvg+tmpSD) && tmp<(tmpAvg+2*tmpSD)] <- tmpAvg+2*tmpSD;
        
        # Initialize current slice properties 
        slice    <- tmp[rBeg:rEnd,cBeg:cEnd];
        sliceAvg <- mean(slice);
        SliceLB  <- tmpAvg + (2 * tmpSD);
        SliceUB  <- tmpAvg + (6 * tmpSD);
        
        # Determine approximate each slice with cylindar
        while( sliceAvg > SliceLB && sliceAvg < SliceUB ) {
          slice    <- tmp[rBeg:rEnd,cBeg:cEnd];
          sliceAvg <- mean(slice);
          
          # Expand cirlce bounduary
          rBeg = rBeg - 1;
          rEnd = rEnd + 1;
          cBeg = cBeg - 1;
          cEnd = cEnd + 1;
        }
        
        # Estimate ventricular radius & volume of this cylindar-ish cross-section
        radEst <- mean(0.5*(rEnd-rBeg)*spacing[1], 0.5*(cEnd-cBeg)*spacing[2]);
        #sliceRad <- cbind(sliceRad, radEst);
        sliceRad[index] <- radEst;
      }
      
      print(sprintf("  => Estimate volume "));
      df0 <- cbind(df0,sliceRad);
      df0_colnames <- c(df0_colnames,sax);
      
      #      sliceVol <- NULL;
      # V = (pi*h/3)*(R1**2+R1*R2+R2**2) = Volume of circular cone frustrum
      #      sliceVol <- cbind(sliceVol, pi*sliceRad*sliceRad*SliceThickness[index]*0.001);    
      #     row.names(sliceVol) <- sax;
      #      if( dfCol != dim(sliceVol)[2] ) {
      #        dfi   <- dfi + 1;
      #        dfCol <- dim(sliceVol)[2];
      #      }
      #      df <- rbind(df,sliceVol);
      
    }
  }
  
  df0[df0==0] <- NA;
  na.approx(df0, 1:length(dcm$img)); 
  colnames(df0) <- df0_colnames;
  rownames(df0) <- rownames(df0, do.NULL=FALSE, prefix="obs.");
  
  
  # Determine the average ejection fraction based on sex and age.  Based on:
  # "Age and gender specific normal values of left ventricular mass, volume and function for gradient echo magnetic resonance imaging: a cross sectional study",
  # By Peter A Cain, Ragnhild Ahl, Erik Hedstrom, Martin Ugander, Ase Allansdotter-Johnsson, Peter Friberg, Hakan Arheden
  # BMC Med Imaging. 2009; 9: 2. Published online 2009 January 21. doi: 10.1186/1471-2342-9-2, PMCID: PMC2657902
  EF_Normal <- 67.0;
  age <- as.integer(substr(PatientsAge,1,nchar(PatientsAge)-1));
  if( PatientsSex=="M" ) {
    if( age >10 && age <21 )        { EF_Normal <- 67.0;  
    } else if( age >20 && age <31 ) { EF_Normal <- 66.0;
    } else if( age >30 && age <51 ) { EF_Normal <- 65.0;
    } else if( age >50 && age <71 ) { EF_Normal <- 64.0;
    } else if( age >70 && age <81 ) { EF_Normal <- 63.0;
    } else if( age >80 )            { EF_Normal <- 62.0; }
  } else if( PatientsSex=="F" ) {
    if( age >10 && age <21 )        { EF_Normal <- 71.0; 
    } else if( age >20 && age <41 ) { EF_Normal <- 70.0;
    } else if( age >40 && age <71 ) { EF_Normal <- 69.0;
    } else if( age >70 && age <81 ) { EF_Normal <- 68.0; 
    } else if( age >80 )            { EF_Normal <- 67.0; } 
  }
  
  # Calculate ejection fraction based on measurements taken from images and 
  # characterize estimate with a diagnosis.
  Diastole <- max(colSums(df));
  Systole  <- min(colSums(df));
  EF <- 100 * (Diastole-Systole) / Diastole;  
  if( EF >= 75 ) {                  Diagnosis <- "Hyperdynamic"; 
  } else if( EF >= 65 && EF < 75) { Diagnosis <- "Normal to Hyperdynamic"; 
  } else if( EF >= 55 && EF < 65) { Diagnosis <- "Normal EF"; 
  } else if( EF >= 45 && EF < 55) { Diagnosis <- "Mildly Abnormal"; 
  } else if( EF >= 35 && EF < 45) { Diagnosis <- "Moderately Abnormal"; 
  } else if( EF >   0 && EF < 35) { Diagnosis <- "Seriously Abnormal"; 
  } else {                          Diagnosis <- "Dead On Arrival"; }
  
  final <- rbind(final,cbind(PatientID,PatientsName,PatientsSex,PatientsAge,Diastole,Systole,EF,Diagnosis))
  print(sprintf("ANALYSIS | Patient #%s", PatientID));  
  print(sprintf("-----------------------------------"));  
  print(sprintf("  Short Axis Stack Size .: %d",    short_axial_stack_size));  
  print(sprintf("  Images Processed ......: %d",    number_of_images));  
  print(sprintf("  Patient Name ..........: %s",    PatientsName));  
  print(sprintf("  Patient Sex ...........: %s",    PatientsSex));  
  print(sprintf("  Patient Age ...........: %s",    PatientsAge));  
  print(sprintf("  Diastole Volume .......: %f ml", Diastole));
  print(sprintf("  Systole Volume ........: %f ml", Systole));
  print(sprintf("  Ejection Fraction"));
  print(sprintf("    Normal ..............: %f %%", EF_Normal));
  print(sprintf("    Calculated ..........: %f %%", EF));
  print(sprintf("    Diagnosis ...........: %s",    Diagnosis));
  print(sprintf("-----------------------------------"));
  print(" ");
  
}
#8-27-34-3-19-10

