library(oro.dicom);
library(oro.nifti);
library(stats);
library(spatstat);
library(mmand);
library(plotrix);
library(seriation);
library(compiler);
library(fields);
library(zoo);
library(xlsx);

############################
## User-Defined Functions ##
############################

#' @name        fn_FitEllipse
#' @author      John Minter
#' @date        2012 
#' @description Calculate ejection fraction (EF) from and characterize 
#'              Least squares fitting of an ellipse to point data 
#' @references  Using the algorithm described in: 
#'                  Radim Halir & Jan Flusser. 1998. 
#'                  Numerically stable direct least squares fitting of ellipses. 
#'                  Proceedings of the 6th International Conference in Central Europe 
#'                  on Computer Graphics and Visualization. WSCG '98, p. 125-132 
#'              Adapted from the original Matlab code by Michael Bedward (2010)
#'              michael.bedward@gmail.com.  Subsequently improved by John Minter (2012)
#'              http://r.789695.n4.nabble.com/Fitting-a-half-ellipse-curve-tp2719037p2720560.html
#' @param  x, y - x and y coordinates of the data points. If a single arg  
#'         is provided it is assumed to be atwo column matrix.
#' @return A list with the following elements:  
#'         coef   - coefficients of the ellipse as described by the general  
#'                  quadratic:  ax^2 + bxy + cy^2 + dx + ey + f = 0  
#'         center - center x and y 
#'         major  - major semi-axis length 
#'         minor  - minor semi-axis length
#'            
fn_FitEllipse <- cmpfun(function(x, y=NULL) {
  
  EPS <- 1.0e-8;
  dat <- xy.coords(x, y); 
  
  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y); 
  D2 <- cbind(dat$x, dat$y, 1); 
  S1 <- t(D1) %*% D1; 
  S2 <- t(D1) %*% D2; 
  S3 <- t(D2) %*% D2; 
  T <- -solve(S3) %*% t(S2); 
  M <- S1 + S2 %*% T; 
  M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2); 
  evec <- eigen(M)$vec; 
  cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2; 
  a1 <- evec[, which(cond > 0)]; 
  f <- c(a1, T %*% a1); 
  names(f) <- letters[1:6];
  
  # calculate the center and lengths of the semi-axes 
  #
  # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2288654/
  # J. R. Minter
  # for the center, linear algebra to the rescue
  # center is the solution to the pair of equations
  # 2ax +  by + d = 0
  # bx  + 2cy + e = 0
  # or
  # | 2a   b |   |x|   |-d|
  # |  b  2c | * |y| = |-e|
  # or
  # A x = b
  # or
  # x = Ainv b
  # or
  # x = solve(A) %*% b
  A <- matrix(c(2*f[1], f[2], f[2], 2*f[3]), nrow=2, ncol=2, byrow=T );
  b <- matrix(c(-f[4], -f[5]), nrow=2, ncol=1, byrow=T);
  soln <- solve(A) %*% b;
  
  b2 <- f[2]^2 / 4;
  
  center <- c(soln[1], soln[2]); 
  names(center) <- c("x", "y"); 
  
  num  <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6]); 
  den1 <- (b2 - f[1]*f[3]); 
  den2 <- sqrt((f[1] - f[3])^2 + 4*b2); 
  den3 <- f[1] + f[3]; 
  
  semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3)) )) 
  
  # calculate the angle of rotation 
  term <- (f[1] - f[3]) / f[2]; 
  angle <- atan(1 / term) / 2; 
  
  return (list(coef=f, center=center, major=max(semi.axes), minor=min(semi.axes), angle=unname(angle))); 
  
});


#' @name         fn_Diagnosis
#' @author       Steven Rankine, Systems Engineer
#'               Kendall-Greene Associates
#' @date         2016-01-28
#' @description  Calculate ejection fraction (EF) from and characterize 
#'               estimate with a diagnosis.
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
#' 
fn_Diagnosis <- cmpfun(function(Diastole, Systole, Age, Sex) {
  
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
  
  # Calculate Ejection Fraction
  Diastole <- round(Diastole, digits=1);
  Systole  <- round(Systole, digits=1);
  if(Diastole==0) {
    EF <- 500.0; 
  } else {
    EF <- 100.0 * (Diastole-Systole) / Diastole; 
  }
  EF <- round(EF, digits=3);
  
  # Determine descriptive diagnosis, based on range that EF falls into
  iEF <- as.integer(EF); # Round EF to the nearesr integer 
  if(iEF == 500)          { Diagnosis <- "Undetermined"; } 
  else if(iEF > 75)       { Diagnosis <- "Hyperdynamic"; } 
  else if(iEF %in% 65:74) { Diagnosis <- "Normal to Hyperdynamic"; } 
  else if(iEF %in% 55:64) { Diagnosis <- "Normal EF"; } 
  else if(iEF %in% 45:54) { Diagnosis <- "Mildly Abnormal"; } 
  else if(iEF %in% 35:44) { Diagnosis <- "Moderately Abnormal"; } 
  else if(iEF %in%  0:34) { Diagnosis <- "Seriously Abnormal"; }
  else                    { Diagnosis <- "Dead On Arrival"; }
  
  return (list(NormalEF=EFNorm, EF=EF, Diagnosis=Diagnosis, Diastole=Diastole, Systole=Systole));
  
});


#' @name        fn_Centroid
#' @author      Steven Rankine, Systems Engineer
#'              Kendall-Greene Associates
#' @date        2016-03-16
#' @description  
#'          
#' @param  img   -  
#' @return   
#' 
fn_Centroid <- cmpfun(function(img) {
  
  mx <- matrix(1:dim(img)[1],nrow=dim(img)[1],ncol=dim(img)[2],byrow=TRUE) *img; 
  my <- matrix(1:dim(img)[2],nrow=dim(img)[1],ncol=dim(img)[2],byrow=FALSE)*img; 
  m  <- sum(img[img!=0]);
  cx <- as.integer(sum(mx[mx!=0])/m);
  cy <- as.integer(sum(my[my!=0])/m);
  return (c(cx,cy));
  
});


#' @name        fn_DCBias
#' @author      Steven Rankine, Systems Engineer
#'              Kendall-Greene Associates
#' @date        2016-03-16
#' @description  
#'          
#' @param  df_in   -  
#' @return   
#' 
fn_DCBias <- cmpfun(function(imgs, Rows, Columns) {
  
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
  
  tmp <- result$roi;
  tmp[result$roi < (mean(result$roi) + 2.5*sd(result$roi))] <- 0;
  result$roi.centroid <- fn_Centroid(tmp);
  return (result);
  
});


#' @name        fn_FrustrumEllipticalVolume
#' @author      Steven Rankine, Systems Engineer
#'              Kendall-Greene Associates
#' @date        2016-02-25
#' @description This function calculate the Volume in milli-liters
#'              of elliptical cone frustrum for each segment. 
#' @param a1 - major axis of upper ellipse of the frustum in mm
#' @param b1 = minor axis of upper ellipse of the frustum in mm
#' @param a2 = major axis of lower ellipse of the frustum in mm
#' @param b2 = minor axis of lower ellipse of the frustum in mm
#' @param h = height of frustum in mm
#' @return  Frustrum elliptical volume in ml
#' 
fn_EllipticalFrustrumVolume <- cmpfun(function( a1, b1, a2, b2, h ) {
  
  area_top <- pi * a1 * b1;
  area_bot <- pi * a2 * b2;
  f_volume <- (h/3)*(area_top + sqrt(area_top*area_bot) + area_bot);
  f_volume <- round(f_volume*0.001, digits=3);
  return (f_volume);
  
});


#' @name        fn_RandomSearch
#' @author      Steven Rankine, Systems Engineer
#'              Kendall-Greene Associates
#' @date        2016-03-16
#' @description  
#'          
#' @param  img     -  
#' @param  ctr     -  
#' @param  rthres  -  
#' @return   
#' 
fn_RandomSearch <- cmpfun(function(img, ctr) {
  
  d <- 64;
  test1 <- ctr[1]; test2 <- ctr[2]; 
  pass  <- 1; pass_max <- as.integer(dim(img)[1]*dim(img)[2]/(d**2));
  while(pass < pass_max){
    n1 <- as.integer(dim(img)[1]/d); s1 <- (ctr[1]-n1):(ctr[1]+n1);
    p1 <- c(rep(2/length(s1),n1),0,rep(0.5/length(s1),n1));
    test1 <- sample(s1,size=1,prob=p1);
    n2 <- as.integer(dim(img)[2]/d); s2 <- (ctr[2]-n2):(ctr[2]+n2);
    p2 <- c(rep(0.5/length(s2),n2),0,rep(2/length(s2),n2));
    test2 <- sample(s2,size=1,prob=p2);
    if(test1 %in% 1:dim(img)[1] & test2 %in% 1:dim(img)[2]){
      if(img[test1,test2]==1 & sum(img[test1-5:test1+5,test2-5:test2+5])>=80) {
        ctr[1] <- test1;
        ctr[2] <- test2;
        break;
      } else {
        if(pass == pass_max-1){
          d <- as.integer(d/2);
          pass <- 1;
          pass_max <- as.integer(dim(img)[1]*dim(img)[2]/(d**2));
        } else {
          pass  <- pass + 1;
        }
      }
    }
  }
  return (ctr);
  
});


#' @name        fn_RemoveNaN
#' @author      Steven Rankine, Systems Engineer
#'              Kendall-Greene Associates
#' @date        2016-03-16
#' @description  
#'          
#' @param  df_in   -  
#' @return   
#' 
fn_RemoveNaN <- cmpfun(function(df_in) {
  
  df_out <- df_in;
  last_row <- dim(df_in)[1];
  last_col <- dim(df_in)[2];
  
  for(i in 1:last_row) {
    if(i==1 | i==last_row) {
      df_out[i,] <- na.approx(df_out[i,], na.rm=FALSE);
    }
    if(TRUE %in% (df_out[i,]<0)) {
      df_out[i,df_out[i,]<0] <- NaN;
    }
  }
  
  if(!(FALSE %in% is.na(df_out[,last_col]))) {
    df_out[1,last_col]                      <- df_out[1,last_col-1];
    df_out[as.integer(last_row/2),last_col] <- df_out[as.integer(last_row/2),last_col-1];
    df_out[last_row,last_col]               <- df_out[last_row,last_col-1];
  }
  
  if(!(FALSE %in% is.na(df_out[,1]))) {
    df_out[1,1]                      <- df_out[1,2];
    df_out[as.integer(last_row/2),1] <- df_out[as.integer(last_row/2),2];
    df_out[last_row,1]               <- df_out[last_row,2];
  }
  
  for(j in 1:last_col) {
    df_out[,j] <- na.approx(df_out[,j], na.rm=FALSE);
  }
  
  df_out <- round(df_out, digits=3);
  return (df_out);
  
});


fn_EllipseArea <- cmpfun(function(a,b) {
  return (pi*a*b);
});
fn_EllipseAreaTriangle <- cmpfun(function(a,b,theta1,theta2) {
  r1 <- sqrt((a*cos(theta1))**2+(b*sin(theta1))**2);
  r2 <- sqrt((a*cos(theta2))**2+(b*sin(theta2))**2);
  return (0.5*(r1*r2)*abs(sin(theta1-theta2)));
});
fn_EllipseAreaSector <- cmpfun(function(a,b,theta1,theta2) {
  f1  <- 0.5*a*b*(theta1-atan((b-a)*sin(2*theta1)/((b+a)-(b-a)*sin(2*theta1))));
  f2  <- 0.5*a*b*(theta2-atan((b-a)*sin(2*theta2)/((b+a)-(b-a)*sin(2*theta2))));
  return (f2[2]-f2[1]);
});
fn_EllipsoidFilled <- cmpfun(function(img, x, y, a, b){
  black <- 0.0; white <- 0.0;
  for(h in 1:dim(img)[1]){
    for(k in 1:dim(img)[2]){
      if((((x-h)**2)/(a**2))+(((y-k)**2)/(b**2)) <= 1){
        if(img[h,k]==0) {
          black <- black + 1.0;
        } else {
          white <- white + 1.0;
        }     
      }
    }
  }
  return (white/(white+black));
});


fn_FindEllipse <- cmpfun(function(img, r_thres=3){ 
  
  vec_thres <- 3;
  for(pass in 1:2){
    
    # Loop over rows of image
    rng <- switch(pass,seq(from=1, to=dim(img)[2], by=3),df$y);  
    for(ym in rng){ 
      ibeg <- 1; vec <- which(img[,ym]==1);
      for(iend in 1:length(vec)){
        if(iend>1){
          x1 <- vec[ibeg]; x2 <- vec[iend-1]; x3 <- vec[iend];
          if((x3-x2)>vec_thres || iend==length(vec)){
            r <- ceiling((x2-x1)*0.5);
            if(r > r_thres) {
              xm <- floor(mean(vec[ibeg:iend-1]));
              z  <- img[xm,ym];
              if(is.null(df)){
                df <- data.frame(x=xm,y=ym,z=z,a=r,b=0,h=0,sl=0,sf=0,sv=0,ct=1);
              } else {
                f <- df$x==xm & df$y==ym;
                if(dim(df[f,])[1]==0){
                  df <- rbind(df, data.frame(x=xm,y=ym,z=z,a=r,b=0,h=0,sl=0,sf=0,sv=0,ct=1));
                } else{
                  df[f,"ct"] <- df[f,"ct"] + 1;
                  df[f,"a"]  <- max(df[f,"a"], r);
                }
              }
            }
            ibeg <- iend;
          }
        } 
      }
    }
    
    # Loop over columns of image
    rng <- switch(pass,seq(from=1, to=dim(img)[1], by=3),df$x);  
    for(xm in rng){
      ibeg <- 1; vec <- which(img[xm,]==1);
      for(iend in 1:length(vec)){
        if(iend>1){
          y1 <- vec[ibeg]; y2 <- vec[iend-1]; y3 <- vec[iend];
          if((y3-y2)>vec_thres || iend==length(vec)){
            r  <- ceiling((y2-y1)*0.5);
            if(r > r_thres){
              ym <- floor(mean(vec[ibeg:iend-1]));
              z  <- img[xm,ym];
              if(is.null(df)){
                df <- data.frame(x=xm,y=ym,z=z,a=0,b=r,h=0,sl=0,sf=0,sv=0,ct=1);
              } else {
                f <- df$y==ym & df$x==xm;
                if(dim(df[f,])[1]==0){
                  df <- rbind(df, data.frame(x=xm,y=ym,z=z,a=0,b=r,h=0,sl=0,sf=0,sv=0,ct=1));
                } else{
                  df[f,"ct"] <- df[f,"ct"] + 1;
                  df[f,"b"]  <- max(df[f,"b"], r);
                }
              }
            }
            ibeg <- iend;
          }
        }
      }
    }
  }
  return (df);
});


fn_Reduce <- cmpfun(function(df_in, dt=7){ 
  
  df <- df_in[df_in$a>0 | df_in$b>0,];
  df <- df[order(df$x,df$y),];  
  row.names(df) <- 1:dim(df)[1];
  delta <- 1;
  while(delta != 0) {
    delta <- dim(df)[1];
    for(i in 2:dim(df)[1]){
      j <- i-1;
      dx <- abs(df$x[i]-df$x[j]); 
      dy <- abs(df$y[i]-df$y[j]);
      if(dx<dt & dy<dt){
        df$x[j:i]  <- as.integer(mean(df$x[j:i], na.rm=T));
        df$y[j:i]  <- as.integer(mean(df$y[j:i], na.rm=T));
        df$a[j:i]  <- as.integer(mean(df$a[j:i], na.rm=T));
        df$b[j:i]  <- as.integer(mean(df$b[j:i], na.rm=T));
        df$h[j:i]  <- as.integer(mean(df$h[j:i], na.rm=T));
        df$ct[j:i] <- as.integer(mean(df$ct[j:i],na.rm=T));
        df$sl[j:i] <- round(mean(df$sl[j:i], na.rm=T), digits=3);
        df$sf[j:i] <- round(mean(df$sf[j:i], na.rm=T), digits=3);
        df$sv[j:i] <- round(mean(df$sv[j:i], na.rm=T), digits=3);
      }
    } 
    df <- unique(df[complete.cases(df),]); 
    row.names(df) <- 1:dim(df)[1];
    delta <- delta - dim(df)[1];
  }
  #df <- df[df$ct>mean(df$ct) & df$a>0 & df$b>0,];
  return (df);
  
}); 


fn_ReduceOverlap <- cmpfun(function(df_in){
  
  df <- df_in[complete.cases(df_in),];
  df <- df[df$a>0 & df$b>0,];
  #   imax <- 0; jmax <- 1;
#   while(imax != jmax){
#     imax <- dim(df)[1];
    for (i in 1:dim(df)[1]){
      x1 <- df[i,"x"]; y1 <- df[i,"y"];
      a1 <- df[i,"a"]; b1 <- df[i,"b"];
      ext1 <- c(x1+a1, y1+b1, x1-a1, y1-b1);
      for (j in 1:dim(df)[1]){
        x2 <- df[j,"x"]; y2 <- df[j,"y"];
        a2 <- df[j,"a"]; b2 <- df[j,"b"];
        f2 <- (a2*b2)/(a1*b1+a2*b2); f1 <- (a1*b1)/(a1*b1+a2*b2);
        ext2 <- c(x2+a2, y2+b2, x2-a2, y2-b2);
        if((i != j) & !(NA %in% ext1) & !(NA %in% ext2)){
          a  <- as.integer((df[i,"a"] *f1)+(df[j,"a"] *f2));
          b  <- as.integer((df[i,"b"] *f1)+(df[j,"b"] *f2));
          h  <- as.integer((df[i,"h"] *f1)+(df[j,"h"] *f2));
          ct <- as.integer((df[i,"ct"]*f1)+(df[j,"ct"]*f2));
          sl <- round((df[i,"sl"]*f1)+(df[j,"sl"]*f2), digits=3);
          sf <- round((df[i,"sf"]*f1)+(df[j,"sf"]*f2), digits=3);
          sv <- round((df[i,"sv"]*f1)+(df[j,"sv"]*f2), digits=3);
          if((ext2[1]<=ext1[1]) & (ext2[2]<=ext1[2]) & (ext1[3]<=ext2[3]) & (ext1[4]<=ext2[4])){
              df[i,"x"]  <- df[j,"x"];
              df[i,"y"]  <- df[j,"y"];
              df[i,"a"]  <- a;  df[j,"a"]  <- a;
              df[i,"b"]  <- b;  df[j,"b"]  <- b;
              df[i,"h"]  <- h;  df[j,"h"]  <- h;
              df[i,"sl"] <- sl; df[j,"sl"] <- sl;
              df[i,"sf"] <- sf; df[j,"sf"] <- sf;
              df[i,"sv"] <- sv; df[j,"sv"] <- sv;
              df[i,"ct"] <- ct; df[j,"ct"] <- ct;
          } else if(ext1[1]<=ext2[1] & ext1[2]<=ext2[2] & ext2[3]<=ext1[3] & ext2[4]<=ext1[4]){
            df[j,"x"]  <- df[i,"x"];
            df[j,"y"]  <- df[i,"y"];
            df[j,"a"]  <- a;  df[i,"a"]  <- a;
            df[j,"b"]  <- b;  df[i,"b"]  <- b;
            df[j,"h"]  <- h;  df[i,"h"]  <- h;
            df[j,"ct"] <- ct; df[i,"ct"] <- ct;
            df[j,"sl"] <- sl; df[i,"sl"] <- sl;
            df[j,"sf"] <- sf; df[i,"sf"] <- sf;
            df[j,"sv"] <- sv; df[i,"sv"] <- sv;
          }
        }
      }
    }
    df <- df[complete.cases(df),];
    df <- df[!duplicated(df[1:6]),]; 
#     jmax <- dim(df)[1];
#   }
  return (df);
  
});


#' @name        fn_RankineCircle
#' @author      Steven Rankine, Systems Engineer
#'              Kendall-Greene Associates
#' @date        2016-02-27
#' @description Feature extraction technique for detecting circles and 
#'              ellipses in a binary image. 
#'          
#' @param  img     -  
#' @param  img_var -  
#' @param  rthres  -  
#' @return Data frame  
#' @export
#' 
fn_RankineCircle <- cmpfun(function(img, img_var, r_thres=3){
  
  r_thres <- 3;

  # 
  df <- NULL;
  df <- fn_FindEllipse(img);
  df <- fn_Reduce(df);
  df <- fn_ReduceOverlap(df);

  # Circle Hough Transform Scoring
  f <- df$z==1 & df$a>0 & df$b>0;
  for(i in which(f)) {
    x <- df[i,"x"]; y <- df[i,"y"];
    a <- df[i,"a"]; b <- df[i,"b"];
    xpts <- NULL;   ypts <- NULL;
    score_var <- img_var[x,y]/max(img_var);
    score_fill<- fn_EllipsoidFilled(img,x,y,a,b);
    for(theta in seq(from=0,to=359,by=10)) {
      theta <- theta * pi / 180;
      xm <- as.integer(a*cos(theta))+x;
      ym <- as.integer(b*sin(theta))+y;
      if(xm %in% 1:dim(img)[1] && ym %in% 1:dim(img)[1]){
        score_loc <- 1.0 - sqrt(((x-xm)/max(df$x))**2 + ((y-ym)/max(df$y))**2);
        if(img[xm,ym]==0) {
          df[i,"h"]  <- df[i,"h"] + 1;              # Accumulator
          df[i,"sl"] <- round(score_loc, digits=3); # 
          df[i,"sv"] <- round(score_var, digits=3); # 
          df[i,"sf"] <- round(score_fill,digits=3); # 
          xpts <- c(xpts,xm);
          ypts <- c(ypts,ym);
        }
      }
    }
#     if(!is.null(xpts) & !is.null(ypts)){
#       ellipParam <- fn_FitEllipse(xpts,ypts);
#       df[i,"x"]   <- as.integer(ellipParam$center[1]); 
#       df[i,"y"]   <- as.integer(ellipParam$center[2]); 
#       df[i,"a"]   <- as.integer(ellipParam$major); 
#       df[i,"b"]   <- as.integer(ellipParam$minor); 
#     }
  }
  return (df[f,]);
  
}); 


##################
## Main Program ##
##################


set.seed(2030);
result_stage3 <- NULL;
actual <- read.csv("../data/train.csv");  # read csv file 
for( StudyID in 234:234 ) { 
  # for( StudyID in sample(1:250,250) ) { 
  
  baseDir <- paste0("../data/train/",StudyID,"/study");

  result_stage2 <- NULL;
  result_stage2_cols <- NULL;
  
  print(baseDir);
  for(sax in list.files(path=baseDir, pattern="sax_*")) {

    print(sprintf("  + %s",sax));
    result_stage1 <- NULL;
    
    # Copy files to a temp directory
    dir.create(paste0(baseDir,"/",sax,"/temp"));
    for( fn in list.files(path=paste0(baseDir,"/",sax)) ) {
      # Ignore zero length files
      if(file.info(paste0(baseDir,"/",sax,"/",fn))$size>0) {
        file.copy(from=paste0(baseDir,"/",sax,"/",fn), to=paste0(baseDir,"/",sax,"/temp/",fn));
      }
    }
    dcm <- readDICOM(paste0(baseDir,"/",sax,"/temp"), verbose=FALSE, recursive=FALSE);
    dcm.info <- dicomTable(dcm$hdr);
    unlink((paste0(baseDir,"/",sax,"/temp")), recursive=TRUE);
    
    # 
    Patient.Name      <- unique(dcm.info["0010-0010-PatientsName"][[1]]);            # Patient's full name.
    Patient.Sex       <- unique(dcm.info["0010-0040-PatientsSex"][[1]]);             # Sex of the named patient. Enumerated Values: M = male F = female O = other
    Patient.Age       <- unique(dcm.info["0010-1010-PatientsAge"][[1]]);             # Age of the Patient.
    Rows              <- as.integer(  dcm.info["0028-0010-Rows"][[1]]);              # Number of rows in the image
    Cols              <- as.integer(  dcm.info["0028-0011-Columns"][[1]]);           # Number of columns in the image
    PixelSpacing      <- as.character(dcm.info["0028-0030-PixelSpacing"][[1]]);      # Physical distance in the patient between the center of each pixel (mm)
    SliceThickness    <- as.integer(  dcm.info["0018-0050-SliceThickness"][[1]]);    # Nominal slice thickness, in mm.
    InstanceNumber    <- as.integer(  dcm.info["0020-0013-InstanceNumber"][[1]]);    # A number that identifies this image. 
    SeriesDescription <- as.character(dcm.info["0008-103E-SeriesDescription"][[1]]); # 
    TriggerTime       <- as.character(dcm.info["0018-1060-TriggerTime"][[1]]);       # 	Time interval measured in msec from the start of the R-wave to the beginning of data taking
    NominalInterval   <- as.character(dcm.info["0018-1062-NominalInterval"][[1]]);   # 
    
    filterSax <- SeriesDescription == "sax" & duplicated(InstanceNumber);#
    if(TRUE %in% filterSax) {
      
      dcm$img        <- dcm$img[filterSax];
      Rows           <- Rows[filterSax];
      Cols           <- Cols[filterSax];
      PixelSpacing   <- PixelSpacing[filterSax];
      InstanceNumber <- InstanceNumber[filterSax];
      SliceThickness <- SliceThickness[filterSax];
      TriggerTime    <- TriggerTime[filterSax];
      NominalInterval<- NominalInterval[filterSax];
      
      XCenter        <- matrix(rep(NA,each=length(dcm$img)));
      YCenter        <- matrix(rep(NA,each=length(dcm$img)));
      MajorRadius    <- matrix(rep(NA,each=length(dcm$img)));
      MinorRadius    <- matrix(rep(NA,each=length(dcm$img)));

      # 
      # print("    Analyzing Files");
      dcm.anal <- fn_DCBias(dcm$img, Rows, Cols);
      for(index in 1:length(dcm$img)) {
        
        spacing <- as.numeric(strsplit(PixelSpacing[index]," ")[[1]]);
        
        # Get the next Image
        tmp <- dcm$img[[index]];
        sd_fac  <- 2.5;
        tmp_roi <- tmp;
        while((length(tmp_roi[tmp_roi==0])/length(tmp_roi)) > 0.05 && sd_fac >= 0.0){
          tmp_roi <- tmp;
          tmp_roi[dcm.anal$roi < (mean(dcm.anal$roi) + sd_fac*sd(dcm.anal$roi))] <- 0;
          sd_fac <- sd_fac - 0.1;
        }
        tmp_crop   <- tmp_roi[which(rowSums(tmp) != 0), which(colSums(tmp) != 0)];
        tmp_lowpas <- gaussianSmooth(tmp_crop, 2);
        tmp_thresh <- threshold(tmp_lowpas, method="kmeans");
        tmp_var    <- dcm.anal$var[which(rowSums(tmp) != 0), which(colSums(tmp) != 0)];
        
        # Find ellipse
        ctr <- fn_Centroid(tmp_lowpas);
        ctr <- fn_RandomSearch(tmp_thresh,ctr);
        a <- mean(c(
          which(rev(tmp_thresh[1:ctr[1],ctr[2]])==0)[1],
          which(tmp_thresh[ctr[1]:dim(tmp_thresh)[1],ctr[2]]==0)[1]
        ), na.rm=TRUE)*spacing[1];
        b <- mean(c(
          which(rev(tmp_thresh[ctr[1],1:ctr[2]]==0))[1],
          which(tmp_thresh[ctr[1],ctr[2]:dim(tmp_thresh)[2]]==0)[1]
        ), na.rm=TRUE)*spacing[2];
        r <- mean(c(a*spacing[1],b*spacing[2]), na.rm=TRUE);

        par(mfrow=c(1, 1));
        image(tmp_thresh, col=gray((0:256)/256));
        abline(h=(0:8)/8.0, col="red", lty="dotted");
        abline(v=(0:8)/8.0, col="red", lty="dotted");
        abline(h=ctr[1]/dim(tmp_thresh)[1],col="green");
        abline(v=ctr[2]/dim(tmp_thresh)[2],col="green");
        draw.ellipse(
          ctr[1]/dim(tmp_thresh)[1],
          ctr[2]/dim(tmp_thresh)[2],
          (a/spacing[1])/dim(tmp_thresh)[1], 
          (b/spacing[2])/dim(tmp_thresh)[2], 
          border="green", lwd=3, 
          main=sprintf("%s | %s | #%d", StudyID, sax, InstanceNumber[index])
        );
        
        # *********** DEBUG ***********
#         
#         df <- fn_RankineCircle(tmp_thresh,tmp_var);
#         f  <- df$a>0 & df$b>0 & df$ct>mean(df$ct) & df$h>mean(df$h) & (df$sl/4 + df$sf/4 + 2*df$sv/4) > mean(df$sl/4 + df$sf/4 + 2*df$sv/4);
#         a  <- mean(df[f,"a"]*spacing[1], na.rm=TRUE);
#         b  <- mean(df[f,"b"]*spacing[2], na.rm=TRUE);
#         r  <- mean(c(a,b), na.rm=TRUE);
# 
#         par(mfrow=c(1, 1));
#         img <- tmp_thresh;
#         image(img, col=gray((0:256)/256));
#         draw.ellipse(
#           (df[f,"x"])/dim(img)[1], 
#           (df[f,"y"])/dim(img)[2], 
#           (df[f,"a"])/dim(img)[1], 
#           (df[f,"b"])/dim(img)[2], 
#           border="green", lwd=3, main=paste0(sax," (",as.character(index),")"));
#         abline(h=(df[f,"y"])/dim(img)[2], v=(df[f,"x"])/dim(img)[1], col="red");
#         #drape.plot(1:Rows[index], 1:Cols[index], tmp_var, border=NA, theta=0, phi=45, main="Spectrogram");
        
        # *********** DEBUG ***********
        
        XCenter[index]     <- ctr[1];
        YCenter[index]     <- ctr[2];
        MajorRadius[index] <- a;
        MinorRadius[index] <- b;
      }
      
      # print("    Record results");
      result_stage1 <- data.frame(
        t = InstanceNumber,
        x = XCenter,
        y = YCenter,
        h = SliceThickness,
        a = MajorRadius,
        b = MinorRadius
      );
      result_stage1 <- result_stage1[order(result_stage1$t),];
      row.names(result_stage1) <- result_stage1$t;
      
      result_stage1 <- cbind(result_stage1,
        alm = predict(lm(a~b+t,data=result_stage1)),
        blm = predict(lm(b~a+t,data=result_stage1)),
        vlm = rep(NA,each=length(result_stage1$a))
      );
      
      h <- result_stage1[1,"h"];
      a1 <- result_stage1[1,"alm"];
      b1 <- result_stage1[1,"blm"];
      result_stage1[1,"vlm"] <- fn_EllipticalFrustrumVolume(a1, b1, a1, b1, h/2);
      for(k in 2:dim(result_stage1)[1]) {
        a2 <- result_stage1[k,"alm"];
        b2 <- result_stage1[k,"blm"];
        result_stage1[k,"vlm"] <- fn_EllipticalFrustrumVolume(a1, b1, a2, b2, h);
        h <- result_stage1[k,"h"];
        a1 <- result_stage1[k,"alm"];
        b1 <- result_stage1[k,"blm"];
      }
    }
    
    if(TRUE %in% filterSax) {
      if(nchar(sax)==5){
        result_stage2_cols <- c(result_stage2_cols, sub("sax_","sax_0",sax));
      } else {
        result_stage2_cols <- c(result_stage2_cols, sax);
      }
      if(is.null(result_stage2)){
        result_stage2 <- data.frame(result_stage1$vlm);
      } else {
        result_stage2 <- cbind(result_stage2,result_stage1$vlm);
      }
      colnames(result_stage2) <- result_stage2_cols;
      result_stage2 <- result_stage2[,order(result_stage2_cols)];   
    }
  }

  # STAGE 2
  result_stage2 <- fn_RemoveNaN(result_stage2);
  result_stage2 <- cbind(result_stage2, Volume=rowSums(result_stage2));
  
  # STAGE 3
  diag.true <- fn_Diagnosis(actual[index,"Diastole"], actual[index,"Systole"], Patient.Age, Patient.Sex);
  diag.calc <- fn_Diagnosis(max(result_stage2[,"Volume"],na.rm=T), min(result_stage2[,"Volume"],na.rm=T), Patient.Age, Patient.Sex);
  nextRow <- cbind(
    Study     = StudyID, 
    Name      = Patient.Name, 
    Sex       = Patient.Sex, 
    Age       = Patient.Age, 
    Ef.Normal = diag.true$NormalEF,
    
    Vd.True   = diag.true$Diastole, 
    Vs.True   = diag.true$Systole,
    Ef.True   = diag.true$EF,
    Diag.True = diag.true$Diagnosis,
    
    Vd.Calc   = diag.calc$Diastole, 
    Vs.Calc   = diag.calc$Systole,
    Ef.Calc   = diag.calc$EF,
    Diag.Calc = diag.calc$Diagnosis
  );
  if(is.null(result_stage3)) {
    result_stage3 <- data.frame(nextRow);
  } else {
    result_stage3 <- rbind(result_stage3,nextRow);
  }
}

write.xlsx(result_stage3, paste0(Sys.Date(),"-1-result_stage3.xlsx"));
