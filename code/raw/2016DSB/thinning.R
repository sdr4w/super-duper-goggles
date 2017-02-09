## http://blog.jakemdrew.com/2013/04/28/image-thinning-using-r/
library(seriation)

absDiff <- function(matrix1,matrix2)
{
  r <- nrow(matrix1)
  c <- ncol(matrix1)
  destMatrix <- matrix1
  for(r in 0:r-1)
  {
    for(c in 0:c-1)
    {
      destMatrix[r,c] <- abs(matrix1[r,c]-matrix1[r,c])
    }
  }
  return(destMatrix)
}

countNonZero <- function(inputMatrix)
{
  return(length(inputMatrix[inputMatrix > 0]))
}

thinningIteration <- function(imageMatrix, iter)
{
  imageInput <- imageMatrix
  r <- nrow(imageInput) - 1
  c <- ncol(imageInput) - 1
  for(i in 2:r)
  {
    for(j in 2:c)
    {
      p2 <- imageInput[i-1, j]
      p3 <- imageInput[i-1, j+1]
      p4 <- imageInput[i, j+1]
      p5 <- imageInput[i+1, j+1]
      p6 <- imageInput[i+1, j]
      p7 <- imageInput[i+1, j-1]
      p8 <- imageInput[i, j-1]
      p9 <- imageInput[i-1, j-1]
      A  <- (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) + 
        (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) + 
        (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
        (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1)
      B  <- p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
      if(iter == 0){
        m1 <- (p2 * p4 * p6)
        m2 <- (p4 * p6 * p8)
      }
      else {
        m1 <- (p2 * p4 * p8)
        m2 <- (p2 * p6 * p8)
      }
      if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
      {
        imageInput[i,j] <- 0
      }
    }
  }
  return(imageInput)
}

thinImage <- function(imageMatrix)
{
  
  im <- imageMatrix
  prev <- im
  repeat {
    im <- thinningIteration(im, 0)
    im <- thinningIteration(im, 1)
    diff <- absDiff(im, prev)
    prev <- im
    if(countNonZero(diff) <= 0)
    {
      break
    }
  } 
  
  return(im)
}



for(ctr in which(img == 1)) {
  df[ctr,"acc"] <- df[ctr,"acc"] + sum(
    which(rev(img[1:df[ctr,"x"], df[ctr,"y"]])==0)[1],
    which(img[df[ctr,"x"]:dim(img)[1], df[ctr,"y"]]==0)[1],
    which(rev(img[df[ctr,"x"], 1:df[ctr,"y"]])==0)[1],
    which(img[df[ctr,"x"], df[ctr,"y"]:dim(img)[2]]==0)[1]
  );
}

for(ym in 1:dim(img)[2]){
  ibeg <- 1;
  ivec <- which(img[,ym]==1);
  for(j in 1:length(ivec)){
    if(j>1){
      if((ivec[j]-ivec[j-1]) != 1){
        xm <- as.integer(mean(ivec[ibeg:j-1]));
        if(length(ibeg:j-1) > 10) {
          df[df$x==xm&df$y==ym,"ct"] <- df[df$x==xm&df$y==ym,"ct"] + 1;
          df[df$x==xm&df$y==ym,"a"]  <- max(df[df$x==xm&df$y==ym,"a"], length(ibeg:j-1));
        }
        ibeg <- j;
      }
    }
  }
}  



for(ctr in which(df[df$a==1,"z"] == 1)) {
  for(theta in seq(from=0,to=2*pi,by=pi/18)) {
    for(r in seq(from=1,to=min(dim(img)),by=4)) {
      xm <- as.integer(r*cos(theta)) + df[ctr,"x"];
      ym <- as.integer(r*sin(theta)) + df[ctr,"y"];
      if(xm %in% 1:dim(img)[1] && ym %in% 1:dim(img)[2]) {
        if(subset(df, x==xm & y==ym)[1,"z"]==0) {
          df[ctr,"a"] <- df[ctr,"a"] + 1;
          #print(sprintf("#%d %3.1f Deg %.2f Radius", ctr, theta*180/pi, r));
        }
      } else {
        break;
      }
    }
  }
}

#n <- 1;
#for(x in 1:dim(img)[1]){
#  for(y in 1:dim(img)[2]){
#    df[n,"x"] <- x;
#    df[n,"y"] <- y;
#    n <- n + 1;
#  }
#}
