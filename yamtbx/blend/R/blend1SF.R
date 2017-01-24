#***********************************************************************************************************
#***********************************************************************************************************
#********* blend1.R
#*********
#********* Copyright (C) 2014 Diamond Light Source & Imperial College London
#*********
#********* Authors: James Foadi & Gwyndaf Evans
#*********
#********* This code is distributed under the BSD license, a copy of which is
#********* included in the root directory of this package.
#************************************************************************************************************
#************************************************************************************************************
# First R module of program BLEND (later to be replaced by C++ code)
# It relies on files produced by blend.cpp.
#
#


# Load libraries and auxiliary files
require(MASS, quietly = TRUE, warn.conflicts = FALSE)  ## For rlm (robust regression)

library("Rcpp")

# Has CCP4 been set up?
ccp4 <- Sys.getenv("CCP4")
if (nchar(ccp4) == 0) stop("You need to set up environment for running ccp4 programs")

ncdist_lib_path <- paste0(ccp4,"/lib/librcpp_ncdist.so")
dyn.load(ncdist_lib_path)
sfdist_lib_path <- paste0(ccp4,"/lib/librcpp_sfdist.so")
dyn.load(sfdist_lib_path)


# Functions


#
# Clean datasets. Basically only accept intensities whose partials add up to
# a fraction between acmin and acmax
cleanData <- function(data,acmin=0.95,acmax=1.05)
{
 data[(data$V5 == "N" & data$V6 < acmin) | (data$V5 == "N" & data$V6 > acmax),5] <- NA
 data <- na.omit(data)

 return(data)
}

#
# Determine if radiation damage detection and cure are to be activated for this data set.
# Minimum number of images to activate procedure is 6; minimum resolution is 5 A.
# Returns TRUE for ok, and FALSE if not ok.
okRadDam <- function(data,minIma=6,minReso=5)
{
 ansR <- TRUE

 # Number of images
 #images <- data$V4
 tmp <- range(data$V4,na.rm=TRUE)
 rima <- tmp[2]-tmp[1]
 #if (length(images) < minIma) ansR <- FALSE
 if (rima < minIma) ansR <- FALSE

 # Resolution range
 resos <- data$V9
 if (1/max(resos) > minReso) ansR <- FALSE

 return(ansR)
}

#
aveIvsRes <- function(nbin,ss,II=NULL,m=min(ss),M=max(ss))
{
 # Resize vectors if needed
 idx <- which(ss >= m & ss <= M)

 # Return NULL if no elements are included inside min and max range
 if (length(idx) == 0)
 {
  return(list(NULL,NULL,NULL,NULL))
 }

 # If there are elements carry on
 ss <- ss[idx]
 if (length(II) > 0) II <- II[idx]

 # Break points
 breaks <- seq(m,M,length=(nbin+1))

 # Histogram
 hst <- hist(ss,breaks=breaks,plot=FALSE)

 # Partition intensity in resolution bins
 if (length(II) > 0)
 {
  aveI <- c()
  sdI  <- c()
  for (i in 1:nbin)
  {
   idx <- which(ss >= breaks[i] & ss < breaks[i+1])
   if (i == nbin) idx <- which(ss >= breaks[i] & ss <= breaks[i+1])
   if (length(idx) > 0)
   {
    aveI <- c(aveI,mean(II[idx],na.rm=TRUE))
    sdI  <- c( sdI,  sd(II[idx],na.rm=TRUE))
   }
   if (length(idx) == 0)
   {
    aveI <- c(aveI,NA)
    sdI <- c(sdI,NA)
   }
  }
 }
 if (length(II) == 0)
 {
  aveI <- NULL
  sdI <- NULL
 }
 
 return(list(hst$counts,hst$mids,aveI,sdI))
}

#
Wplots <- function(data,nbin,nwin=1,m=NA,M=NA)
{
 # Analysis of Wilson's plot shape evolution in time (batch dependency)
 # Batch column is number 4
 # Full/Partial flag is column 5
 # Intensity column is number 7
 # sigI column is number 8
 # Resolution column is number 9

 # Extract batches
 batches <- data[,4]

 # Range of batches across which to investigate plot evolution
 rbat <- sort(unique(batches))
 if (length(rbat) < (2*nwin+1)) stop("Number of batches too small for meaningful modeling")

 # Just to extract range (yrange will be worked out in the loop)
 ltmp <- aveIvsRes(nbin,data[,9])
 xrange <- range(ltmp[[2]])
 if (!is.na(m)) xrange[1] <- m
 if (!is.na(M)) xrange[2] <- M

 # Cycle through increasing number of batches (from early to late time during data collection)
 for (i in (nwin+1):(length(rbat)-nwin))
 {
  # Select resolutions and intensities (standardized or not) within batch range
  ileft <- rbat[(i-nwin)]
  iright <- rbat[(i+nwin)]
  s <- data[data[,4]  >= ileft & data[,4] <= iright,9]
  nI <- data[data[,4] >= ileft & data[,4] <= iright,7]
  
  # Create elements of Wilson's plots
  ltmp <- aveIvsRes(nbin,s,nI,m=xrange[1],M=xrange[2])
  if (i == (nwin+1))
  {
   ss <- matrix(ltmp[[2]],nrow=1,byrow=TRUE)
   II <- matrix(ltmp[[3]],nrow=1,byrow=TRUE)
   sDD <- matrix(ltmp[[4]],nrow=1,byrow=TRUE)
  }
  if (i != (nwin+1))
  {
   ss <- rbind(ss,ltmp[[2]])
   II <- rbind(II,ltmp[[3]])
   sDD <- rbind(sDD,ltmp[[4]])
  }
 }

 return(list(ss=ss,II=II,sDD=sDD))
}

#
mrad <- function(wplots,addbit=0.1)
{
 # Find out when plots start to be empty
 ikey <- checkEmptyPlots(wplots)

 # Resolutions to be added to final list
 ss <- wplots[[1]][1,ikey]

 # Models

 # Exponential
 I0 <- c()
 B <- c()
 RR <- c()
 for (i in ikey)
 {
  y0 <- wplots[[2]][,i]
  minY <- min(y0,na.rm=TRUE)
  if (minY > 0) y <- y0
  if (minY <= 0) y <- y0-minY+addbit
  #y[y == 0] <- NA                # To avoid NaN and -Inf in logarithm of zero
  x <- 1:length(y)
  model <- lm(log(y)~x,na.action=na.omit)
  cff <- unname(coef(model))
  if (minY > 0) I0 <- c(I0,exp(cff[1]))
  if (minY < 0) I0 <- c(I0,exp(cff[1])+minY-addbit)
  B <- c(B,-cff[2])
  RR <- c(RR,summary(model)[[8]][1])
 }
 lista_exponential <- list(s=ss,I0=I0,B=B,RR=RR)

 return(lista_exponential)
}

#
# As many bins of each Wilson plot can be including exclusively NA's,
# this function extracts bin number for all bins with some data
checkEmptyPlots <- function(wplots)
{
 nbin <- dim(wplots$II)[2]
 ikey <- c()
 for (i in 1:nbin)
 {
  nna <- sum(is.na(wplots[[2]][,i]))
  if (nna < dim(wplots$II)[1]) ikey <- c(ikey,i)
  #if (nna == 0) ikey <- c(ikey,i)
 }

 return(ikey)
}

#
# New procedure to find out radiation damage and cut off images
newraddam <- function(wplots)
{
 # Radiation damage routines
 tmp <- mrad(wplots)
 model <- rlm(tmp$B ~ tmp$s,maxit=40)
 coefs <- summary(model)$coefficients[2,1:2]
 attributes(coefs) <- NULL

 ans <- NULL
 if (!is.na(coefs[1]))
 {
  if (coefs[1] < 0) ans <- FALSE
 }
 if (is.na(coefs[1])) ans <- FALSE
 if (!is.na(coefs[2]))
 {
  if (!is.na(coefs[1]))
  {
   if (coefs[2] > 0)
   {
    if (coefs[1]-coefs[2] <= 0) ans <- FALSE
    if (coefs[1]-coefs[2] > 0) ans <- TRUE
   }
  }
 }
 if (is.na(coefs[2])) ans <- FALSE

 # If ans is FALSE return infinity
 if (ans)
 {
  coefs <- model$coefficients
  attributes(coefs) <- NULL
 }
 if (!ans)
 {
  coefs <- c(-Inf,Inf)
 }

 return(coefs)
}

#
# Function returning cutting curve, i.e. curve returning from which image
# one should cut off at each resolution
tfit <- function(s,coefs,fdrop=0.75)
{
 Be <- coefs[2]*s+coefs[1]
 if (Be > 0) te <- -log(fdrop)/Be
 if (Be <= 0) te <- NA

 return(te)
}

# To extract a list with datasets to be merged, as many as the nodes of a dendrogram (htree)
# This list can be subsequently used in "merge_datasets" (in a loop) to extract statistics
# for each node of the dendrogram
treeToList <- function(htree,lres,hres)
{
 # Function to list all groups (and resolution ranges) corresponding to each node in a dendrogram.
 # "htree" is the cluster object (from hclust).
 # lres and hres are named vectors. Names must coincide with labels in the cluster object.


 # Check dataset numbers in tree agree with lres and hres length
 nres <- length(lres)
 if (length(hres) != nres) stop("Bad input to function treeToList: lres and hres have differing lengths")
 tmp <- as.vector(htree$merge)
 tmp <- -tmp[tmp < 0]
 if (length(tmp) != nres) stop("Dendrogram deals with a number of datasets differing from length of lres or hres")

 groups <- list()
 groups_resos <- list()
 itree <- dim(htree$merge)[1]
 for (i in 1:itree)
 {
  il <- htree$merge[i,1]
  ir <- htree$merge[i,2]
  if (il < 0) lab1 <- htree$labels[-il]
  if (ir < 0) lab2 <- htree$labels[-ir]
  if (il > 0) lab1 <- groups[[il]]
  if (ir > 0) lab2 <- groups[[ir]]
  lab <- c(lab1,lab2)
  idx <- which(names(lres) %in% lab)
  Lres <- max(lres[idx])
  idx <- which(names(hres) %in% lab)
  #Hres <- min(hres[idx])
  Hres <- max(hres[idx])
  lab <- as.integer(lab)
  groups <- c(groups,list(lab))
  groups_resos <- c(groups_resos,list(c(Lres,Hres)))
 }

 return(list(groups,groups_resos))
}

#
# Interpolate plot resolution vs signal-to-noise ratio with a 10-degrees
# polynomial and get highest resolution which provides given signal-to-noise
bestReso <- function(xr,yr,snr=1.5)
{
 # Flag to warn user if intensities are weak
 weakFLAG <- FALSE

 # Modeling signal-to-noise
 model <- lm(yr~poly(xr,degree=10))
 xp <- seq(range(xr)[1],range(xr)[2],length=10000)
 yp <- predict(model,list(xr=xp))

 # If data are too weak, raise flag and return full reso
 dus <- yp[yp < snr]
 if (length(dus)/10000 > 0.5)
 {
  weakFLAG <- TRUE
  return(c(weakFLAG,1/xp[10000]))
 }

 # Move from highest to lowest resolution until snr is reached
 ip <- 10000
 while (yp[ip] < snr) ip <- ip-1

 return(c(weakFLAG,1/xp[ip]))
}

#
# Reduce data dimensionality and standardize for cell parameters
nparCell <- function(data,cn,cumprob)
{
 # Normalize data
 model <- prcomp(data,scale=TRUE)
 smod <- summary(model)

 # Choose minimal number of variables that give enough statistical variability
 if (length(model$x[1,]) == 1) npar <- model$x
 if (length(model$x[1,]) > 1)
 {
  idx <-  which(smod$importance[3,] > cumprob)
  npar <- model$x[,1:idx[1]]
 }
 npar <- as.matrix(npar)
 rownames(npar) <- cn

 return(npar)
}

#
# f is a vector containing m sampled values of a reference function f.
# g is a vector containing m sampled vaues of a function g to be scaled close to f.
# This function returns the scaling constant k: new g = k X g
scale_two_functions <- function(f,g)
{
 m <- length(f)
 if (length(g) != m) stop("The two sampled functions need to have the same number of sampled values")

 # k is found as the value making sum(kg[i]-f[i])^2 a minimum
 fsx <- which(is.na(f))
 gsx <- which(is.na(g))
 f[gsx] <- NA
 g[fsx] <- NA
 f <- na.omit(f)
 g <- na.omit(g)
 k=sum(f*g)/sum(g^2)

 return(k)
}

#
# Prototype to produce "npar" type matrix on which cluster analysis can be applied
nparWilson <- function(listBplots,cn,nref=1,cumprob=0.95)
{
 # Raise all plots to make them positive
 listBplots <- raiseWplots(listBplots)

 # Reference structure
 wref <- listBplots[[nref]]

 # Scale plots to line them up with reference structure
 for (i in 1:length(listBplots))
 {
  tmp <- wScale(wref[[3]],listBplots[[i]][[3]])
  listBplots[[i]][[3]] <- tmp
 }

 # Build data frame
 wpar <- matrix(listBplots[[1]][[3]],nrow=1,byrow=TRUE)
 for (i in 2:length(listBplots)) wpar <- rbind(wpar,matrix(listBplots[[i]][[3]],nrow=1,byrow=TRUE))

 # Replace Na's if they (alas!) appear in wpar
 wpar <- interpolateNA(listBplots[[1]][[2]],wpar)

 # Principal Component Analysis to reduce dimensionality
 model <- prcomp(wpar,scale=TRUE)
 smod <- summary(model)
 if (length(model$x[1,]) == 1) npar <- model$x
 if (length(model$x[1,]) > 1)
 {
  idx <-  which(smod$importance[3,] > cumprob)
  npar <- model$x[,1:idx[1]]
 }
 npar <- as.matrix(npar)
 rownames(npar) <- cn

 return(npar)
}

#
# To raise plots in case parts of them are negative (used when argument of a logarithm)
raiseWplots <- function(listWplots,addbit=100)
{
 # Find plots minimum
 mWP <- min(listWplots[[1]][[3]],na.rm=TRUE)
 for (i in 2:length(listWplots))   if (min(listWplots[[i]][[3]],na.rm=TRUE) < mWP) mWP <- min(listWplots[[i]][[3]],na.rm=TRUE)

 # Now raise all plots (if mWP is negative)
 if (!is.na(mWP))
 {
  if (mWP < 0) for (i in 1:length(listWplots)) listWplots[[i]][[3]] <- listWplots[[i]][[3]]-mWP+addbit
 }

 return(listWplots)
}

#
# Scale a function g (g is an array) in relation to f (array with same number of points as g)
# using exponential regression. Both g and f are supposed to be positive, i.e. none of the sampled
# points can be zero or negative.
wScale <- function(f,g)
{
 # x and y for regression
 x <- 1:length(f)
 y <- log(g/f)

 # Modeling
 model <- lm(y~x)

 # Coefficients
 coefs <- model$coefficients
 attributes(coefs) <- NULL

 # Scale g
 g <- g*exp(-coefs[2]*x)/exp(coefs[1])

 return(g)
}

#
# To re-arrange file names as changed by file forR_raddam.dat so that
# thay have the same order as file NEW_list_of_files.dat
reArrangeFilenames <- function(filenames)
{
 idx <- c()
 for (nome in filenames)
 {
  tmp <- strsplit(nome,"_")
  tmp <- strsplit(tmp[[1]][3],".",fixed=TRUE)
  idx <- c(idx,as.integer(tmp[[1]][1]))
 }

 # Build data frame and next order
 idxfile <- data.frame(idx=idx,file=I(filenames))
 newidxfile <- idxfile[order(idxfile$idx),]

 return(newidxfile$file)
}

# Functions to calculate cell variation in angstroms
.dcl <- function(a,b,c,aa,bb,gg)
{
# Input: cell parameters. Output: values useful to all crystallographic
# calculations. These are:
# 1) sa = sin(alpha), sb = sin(beta), sc = sin(gamma)
# 2) ca = cos(alpha), cb = cos(beta). cc = cos(gamma)
# 3) sides of reciprocal cell: ar = a*, br = b*, cr = c*
# 4) sines of angles of reciprocal cell: sar = sin(alpha*), sbr = sin(beta*), scr = sin(gamma*)
# 5) cosines of angles of reciprocal cell: car = cos(alpha*), cbr = cos(beta*), ccr = cos(gamma*)
# 6) Volume of unit cell: V
 aa <- aa*pi/180
 bb <- bb*pi/180
 gg <- gg*pi/180
 sa <- sin(aa)
 sb <- sin(bb)
 sc <- sin(gg)
 ca <- cos(aa)
 cb <- cos(bb)
 cc <- cos(gg)

 # To avoid NaN generated by rounding off errors, use for cell-derived quantities formulas
 # derived previously by computationa crystallographers
 sang <- 0.5*(aa+bb+gg)
 V2 <- sqrt(sin(sang-aa)*sin(sang-bb)*sin(sang-gg)*sin(sang))
 V <- 2*a*b*c*V2
 ar <- b*c*sa/V
 br <- a*c*sb/V
 cr <- a*b*sc/V
 car <- (cb*cc-ca)/(sb*sc)
 cbr <- (ca*cc-cb)/(sa*sc)
 ccr <- (ca*cb-cc)/(sa*sb)
 sar <- sqrt(1-car*car)
 sbr <- sqrt(1-cbr*cbr)
 scr <- sqrt(1-ccr*ccr)
 l <- c(sa,sb,sc,ca,cb,cc,ar,br,cr,sar,sbr,scr,car,cbr,ccr,V)
 names(l) <- c("SIN_ALPHA","SIN_BETA","SIN_GAMMA","COS_ALPHA","COS_BETA","COS_GAMMA","A*","B*","C*","SIN_ALPHA*","SIN_BETA*","SIN_GAMMA*",
               "COS_ALPHA*","COS_BETA*","COS_GAMMA*","V")

 return(l)
}

# Given cell parameters, return the inverse of cell orthogonalization matrix (First choice in Giacovazzo's book)
.triclinic_to_orthogonal_01 <- function(a,b,c,aa,bb,cc)
{
 lp <- .dcl(a,b,c,aa,bb,cc)
 m_1 <- matrix(c(a,0,0,b*lp[6],b*lp[3],0,c*lp[5],-c*lp[2]*lp[13],1/lp[9]),nrow=3,ncol=3,byrow=TRUE)

 return(m_1)
}

# Given cell parameters, return the inverse of cell orthogonalization matrix (MOSFLM choice, second choice in Giacovazzo's book)
.triclinic_to_orthogonal_02 <- function(a,b,c,aa,bb,cc)
{
 lp <- .dcl(a,b,c,aa,bb,cc)
 m_1 <- matrix(c(1/lp[7],-lp[15]/(lp[7]*lp[12]),a*lp[5],0,1/(lp[8]*lp[12]),b*lp[4],0,0,c),nrow=3,ncol=3,byrow=TRUE)

 return(m_1)
}

.dist_oblend <- function(cent1,a1,b1,c1,alpha1,beta1,gamma1,cent2,a2,b2,c2,alpha2,beta2,gamma2)
{
    M1 <- .triclinic_to_orthogonal_01(a1,b1,c1,alpha1,beta1,gamma1)
    M2 <- .triclinic_to_orthogonal_01(a2,b2,c2,alpha2,beta2,gamma2)
    v1 <- M1%*%matrix(c(1,1,1),ncol=1)
    v2 <- M2%*%matrix(c(1,1,1),ncol=1)
    dd <- sqrt((v1[1,1]-v2[1,1])^2+(v1[2,1]-v2[2,1])^2+(v1[3,1]-v2[3,1])^2)
    return(dd)
}

.dist_nblend <- function(cent1,a1,b1,c1,alpha1,beta1,gamma1,cent2,a2,b2,c2,alpha2,beta2,gamma2)
{
    dd <- as.double(.Call("rcpp_ncdist", "P", a1,b1,c1,alpha1,beta1,gamma1,"P",a2,b2,c2,alpha2,beta2,gamma2))
    return(dd)
}

.dist_sfblend <- function(mtz_file1,mtz_file2)
{
    dd <- as.double(.Call("rcpp_sfdist", toString(mtz_file1), toString(mtz_file2)))
    return(dd)
}

# Compute matrix of cross-cells distances
evaluateMaxChange <- function(maindf)
{
 n <- length(maindf[,2])
 dMat <- matrix(nrow=n,ncol=n)
 for (i in 1:n)
 {
  for (j in 1:n)
  {
    dMat[i,j] <- .dist_nblend(maindf[i,8],maindf[i,2],maindf[i,3],maindf[i,4],maindf[i,5],maindf[i,6],maindf[i,7],maindf[j,8],maindf[j,2],maindf[j,3],maindf[j,4],maindf[j,5],maindf[j,6],maindf[j,7])
  }
 }

 return(dMat)
}

+# Compute matrix of cross-dataset distances based on mtz structure factor files
+evaluateSFChange <- function(maindf)
{
    filenames <- read.table("./NEW_list_of_files.dat",as.is=c(1))
    for (ii in 1:length(filenames[,1]))
    {
        stmp <- normalizePath(filenames$V1[ii])
        filenames$V1[ii] <- stmp
    }

    n <- length(maindf[,2])
    dMat <- matrix(nrow=n,ncol=n)
    for (i in 1:n)
    {
        for (j in i:n)
        {
            dMat[i,j] <- .dist_sfblend(filenames$V1[i],filenames$V1[j])
            if (j > i) dMat[j,i] <- dMat[i,j]
        }
    }
    
    return(dMat)
}

# Functions to evaluate max change in linear dimensions for cell parameters

# Diagonals along 3 independent faces of unit cell
faceDiagonals <- function(a,b,c,aa,bb,cc)
{
 dab <- sqrt(a^2+b^2-2*a*b*cos(pi-cc*pi/180))
 dac <- sqrt(a^2+c^2-2*a*c*cos(pi-bb*pi/180))
 dbc <- sqrt(b^2+c^2-2*b*c*cos(pi-aa*pi/180))

 return(c(dab,dac,dbc))
}

distRatio <- function(v)
{
 # Length of vector
 n <- length(v)

 # Create n X n matrix
 m <- matrix(ncol=n,nrow=n)

 # Double loop to fill the matrix
 for (i in 1:n)
 {
  for (j in 1:n)
  {
   dd <- abs(v[i]-v[j])
   mv <- min(v[i],v[j])
   m[i,j] <- (dd/mv)*100
  }
 }

 return(m)
}

adistRatio <- function(v)
{
 # Length of vector
 n <- length(v)

 # Create n X n matrix
 m <- matrix(ncol=n,nrow=n)

 # Double loop to fill the matrix
 for (i in 1:n)
 {
  for (j in 1:n)
  {
   m[i,j] <- abs(v[i]-v[j])
  }
 }

 return(m)
}

# Get matrix indices
get_indices <- function(N,n)
{
 # N is the serial number and n the size of n X n matrix
 i <- N%%n
 if (i == 0) i <- n
 j <- N%/%n+1
 if (i == n) j <- j-1
 if (j > n) j <- n

 return(c(i,j))
}

maxRatio <- function(macropar,idx)
{
 # Cell parameters
 cpar <- macropar[idx,2:7]

 # Number of datasets
 n <- length(cpar[,1])
 
 # 3 diagonal lengths for all datasets
 dab <- c()
 dac <- c()
 dbc <- c()
 for (i in 1:n)
 {
  tmp <- faceDiagonals(cpar[i,1],cpar[i,2],cpar[i,3],cpar[i,4],cpar[i,5],cpar[i,6])
  dab <- c(dab,tmp[1])
  dac <- c(dac,tmp[2])
  dbc <- c(dbc,tmp[3])
 }

 # Calculate maxRatio matrix for the 3 diagonals vectors and extract max value for each matrix
 #mab <- max(distRatio(dab))
 mab <- max(adistRatio(dab))
 #iab <- which(distRatio(dab) == max(distRatio(dab)))
 iab <- which(adistRatio(dab) == max(adistRatio(dab)))
 ijab <- get_indices(iab[1],n)
 #sab <- adistRatio(dab)[iab]
 sab <- distRatio(dab)[iab]
 #mac <- max(distRatio(dac))
 mac <- max(adistRatio(dac))
 #iac <- which(distRatio(dac) == max(distRatio(dac)))
 iac <- which(adistRatio(dac) == max(adistRatio(dac)))
 ijac <- get_indices(iac[1],n)
 #sac <- adistRatio(dac)[iac]
 sac <- distRatio(dac)[iac]
 #mbc <- max(distRatio(dbc))
 mbc <- max(adistRatio(dbc))
 #ibc <- which(distRatio(dbc) == max(distRatio(dbc)))
 ibc <- which(adistRatio(dbc) == max(adistRatio(dbc)))
 ijbc <- get_indices(ibc[1],n)
 #sbc <- adistRatio(dbc)[ibc]
 sbc <- distRatio(dbc)[ibc]
 vv <- c(sab[1],sac[1],sbc[1])
 ij <- matrix(c(ijab,ijac,ijbc),ncol=3)
 imall <- which(c(mab,mac,mbc) == max(c(mab,mac,mbc)))[1]
 Mpar <- vv[imall]
 ijs <- ij[,imall]
 cns <- macropar[idx[ijs],"cn"]

 #return(c(max(mab,mac,mbc),Mpar,cns[1],cns[2]))
 return(c(Mpar,max(mab,mac,mbc),cns[1],cns[2]))
}

find_nodes_coords <- function(clst, clns, cn)
{
 # Number of objects
 nobj <- length(clst$merge[,1]) + 1

 # x coordinates for individual objects
 xobj <- clns[[nobj - 1]]

 # Go through nodes formation
 x <- numeric(nobj - 1)
 y <- numeric(nobj - 1)
 for (i in 1:length(clst$merge[,1]))
 {
  ele1 <- clst$merge[i,1]
  ele2 <- clst$merge[i,2]
  if (ele1 < 0 & ele2 < 0)
  {
   #idx1 <- match(abs(ele1), xobj)
   #idx2 <- match(abs(ele2), xobj)
   idx1 <- match(cn[abs(ele1)], xobj)
   idx2 <- match(cn[abs(ele2)], xobj)
   x[i] <- 0.5 * (idx1 + idx2)
   y[i] <- clst$height[i]
  }
  if (ele1 < 0 & ele2 > 0)
  {
   #idx1 <- match(abs(ele1), xobj)
   idx1 <- match(cn[abs(ele1)], xobj)
   idx2 <- x[ele2]
   x[i] <- 0.5 * (idx1 + idx2)
   y[i] <- clst$height[i]
  }
  if (ele1 > 0 & ele2 < 0)
  {
   idx1 <- x[ele1]
   #idx2 <- match(abs(ele2), xobj)
   idx2 <- match(cn[abs(ele2)], xobj)
   x[i] <- 0.5 * (idx1 + idx2)
   y[i] <- clst$height[i]
  }
  if (ele1 > 0 & ele2 > 0)
  {
   idx1 <- x[ele1]
   idx2 <- x[ele2]
   x[i] <- 0.5 * (idx1 + idx2)
   y[i] <- clst$height[i]
  }
 }

 return(list(x=x,y=y))
}

# The following function, at present, uses a very crude and badly-thought interpolation.
# It needs more investigations and proper fixing. It is used wherever prcomp is called
interpolateNA <- function(x,M)
{
 # Extract rows and columns
 tmp <- dim(M)
 m <- tmp[1]
 n <- tmp[2]

 # Find rows where there are NA's
 na_rows <- c()
 for (i in 1:m)
 {
  jfl <- 0
  for (j in 1:n)
  {
   if (is.na(M[i,j])) jfl <- 1
  }
  if (jfl == 1) na_rows <- c(na_rows,i)
 }

 # Return unchanged M if no NA's found
 if (length(na_rows) == 0) return(M)

 # Mend each "NA's-containing" row with an interpolated value
 for (i in na_rows)
 {
  # Find out positions of NA and their complement
  idx <- which(is.na(M[i,]))
  jdx <- which(!is.na(M[i,]))

  # 10-degree polynomial interpolation (normally degree 10; lower if less points available)
  xr <- x[jdx]
  yr <- M[i,jdx]
  dgree <- min(length(yr),11)-1
  #model <- lm(yr~poly(xr,degree=10))
  model <- lm(yr~poly(xr,degree=dgree))
  xp <- x
  yp <- predict(model,list(xr=xp))

  # If interpolated NA values are negative, replace with 0.1
  for (j in idx) if (yp[j] < 0) yp[j] <- 0.1

  # Replacing NA values in original array
  M[i,idx] <- yp[idx]
 }

 return(M)
}
####################################################################
####################################################################
######################### Main program #############################
####################################################################
####################################################################

# To avoid warning messages set warn to a negative value 
options(warn = -1)

# Fix global parameters here (default values)
nbin <- 20                    # Number of bins for dynamic and overall Wilson plots
fdrop <- 0.75                 # Intensity decay fraction for datasets affected by radiation damage (keyword RADFRAC)
isigi <- 1.5                  # Signal-to-noise ratio wanted to estimate highest resolution cut
cparweight <- 1.0             # Weight to be given to cluster analysis with cell parameters as descriptors.
                              # Wilson plots descriptors get automatically weight 1-cparweight

#######################################################################################################
#######################################################################################################
######################################## KEYWORDS FROM FILE ###########################################
#######################################################################################################
#######################################################################################################
# Check whether there's a ready keyword file and use it if found
if (file.exists("BLEND_KEYWORDS.dat"))
{
 contents <- scan("BLEND_KEYWORDS.dat", what="character",sep="\n",quiet=TRUE)

 # Find section "BLEND KEYWORDS"
 idxs <- grep("BLEND KEYWORDS",contents,fixed=TRUE)
 idxe <- grep("POINTLESS KEYWORDS",contents,fixed=TRUE)

 # Only extract values if something is included in "BLEND KEYWORDS" section
 if ((idxe-idxs) > 1)
 {
  idxs <- idxs+1
  idxe <- idxe-1

  # Turn values into data frame
  slm <- strsplit(contents[idxs:idxe]," ")
  tmp <- data.frame()
  for (slt in slm)
  {
   jdx <- which(nchar(slt) != 0)
   #if (length(jdx) != 2) stop("Wrongly formatted BLEND_KEYWORDS.dat file")
   tmp <- rbind(tmp,data.frame(I(slt[jdx[1]]),I(slt[jdx[2]])))             # "I()" to avoid data frame turning characters into factors
  }

  for (i in 1:length(tmp[,1]))
  {
   if (as.character(tmp[i,1]) == "NBIN" | as.character(tmp[i,1]) == "nbin") nbin <- as.numeric(tmp[i,2])
   if (as.character(tmp[i,1]) == "RADFRAC" | as.character(tmp[i,1]) == "radfrac" |
       as.character(tmp[i,1]) == "RADF"    | as.character(tmp[i,1]) == "radf") fdrop <- as.numeric(tmp[i,2])
   if (as.character(tmp[i,1]) == "ISIGI" | as.character(tmp[i,1]) == "isigi" |
       as.character(tmp[i,1]) == "ISIG"  | as.character(tmp[i,1]) == "isig") isigi <- as.numeric(tmp[i,2])
   if (as.character(tmp[i,1]) == "CPARWT" | as.character(tmp[i,1]) == "cparwt" |
       as.character(tmp[i,1]) == "CPAR"   | as.character(tmp[i,1]) == "cpar") cparweight <- as.numeric(tmp[i,2])
  }

  # Turn indices back to their initial value for following section
  idxs <- idxs-1
  idxe <- idxe+1
 }

 # Not all keywords can be included in existing "BLEND_KEYWORDS.dat". Complete as needed.
 cat("BLEND KEYWORDS\n",file="BLEND_KEYWORDS.dat")
 if ((idxe-idxs) > 1)
 {
  idxs <- idxs+1
  idxe <- idxe-1
  for (i in idxs:idxe) cat(paste(contents[i],"\n",sep=""),file="BLEND_KEYWORDS.dat",append=TRUE)
  nomi <- as.character(tmp[,1])
  if (!("NBIN" %in% nomi) & !("nbin" %in% nomi))
  {
   linea <- sprintf("NBIN      %d\n",as.integer(nbin))
   cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  }
  if (!("RADFRAC" %in% nomi) & !("radfrac" %in% nomi) & !("RADF" %in% nomi) & !("radf" %in% nomi))
  {
   linea <- sprintf("RADFRAC   %5.3f\n",as.numeric(fdrop))
   cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  }
  if (!("ISIGI" %in% nomi) & !("isigi" %in% nomi) & !("ISIG" %in% nomi) & !("isig" %in% nomi))
  {
   linea <- sprintf("ISIGI     %5.3f\n",as.numeric(isigi))
   cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  }
  if (!("CPARWT" %in% nomi) & !("cparwt" %in% nomi) & !("CPAR" %in% nomi) & !("cpar" %in% nomi))
  {
   linea <- sprintf("CPARWT    %5.3f\n",as.numeric(cparweight))
   cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  }

  # Turn indices back to their initial value for following section
  idxs <- idxs-1
  idxe <- idxe+1
 }

 # If nothing was included in BLEND KEYWORDS section add everything
 if ((idxe-idxs) == 1)
 {
  linea <- sprintf("NBIN      %d\n",as.integer(nbin))
  cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  linea <- sprintf("RADFRAC   %5.3f\n",as.numeric(fdrop))
  cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  linea <- sprintf("ISIGI     %5.3f\n",as.numeric(isigi))
  cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
  linea <- sprintf("CPARWT    %5.3f\n",as.numeric(cparweight))
  cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
 }

 # Append rest of keywords
 for (i in idxe:(length(contents)-1)) cat(paste(contents[i],"\n",sep=""),file="BLEND_KEYWORDS.dat",append=TRUE)
 cat(contents[length(contents)],file="BLEND_KEYWORDS.dat",append=TRUE)
}

# If keyword file is not found write one with default words
if (!file.exists("BLEND_KEYWORDS.dat"))
{
 cat("BLEND KEYWORDS\n",file="BLEND_KEYWORDS.dat")
 linea <- sprintf("NBIN      %d\n",as.integer(nbin))
 cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
 linea <- sprintf("RADFRAC   %5.3f\n",as.numeric(fdrop))
 cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
 linea <- sprintf("ISIGI     %5.3f\n",as.numeric(isigi))
 cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)
 linea <- sprintf("CPARWT    %5.3f\n",as.numeric(cparweight))
 cat(linea,file="BLEND_KEYWORDS.dat",append=TRUE)

 # Complete with keywords sections for POINTLESS and AIMLESS
 cat("POINTLESS KEYWORDS\n",file="BLEND_KEYWORDS.dat",append=TRUE)
 cat("AIMLESS KEYWORDS\n",file="BLEND_KEYWORDS.dat",append=TRUE)
}

#######################################################################################################
#######################################################################################################
################################### RADIATION DAMAGE ANALYSIS #########################################
#######################################################################################################
#######################################################################################################

# Load input file names
data <- read.table(file="forR_raddam.dat")
nfiles <- length(data[,1])
filenames <- as.character(data$V1)

# Use ordering of file NEW_list_of_files.dat, rather than order from forR_raddam.dat
filenames <- reArrangeFilenames(filenames)

# Loop to get rid of radiation damaged observations
# and to compute overall, relative B factor
sMax <- c()
sMin <- c()
listData <- list()
batches <- c()
Ibatches <- c()
Fbatches <- c()
cat("Radiation damage analysis ...\n")
tempWplots <- list()
for (crystal in filenames)
{
 cat(sprintf("Dealing with file %s.\n",crystal))
 cat(sprintf("Considering purging data of radiation damaged elements."))
 data <- read.table(crystal)
 data <- cleanData(data)

 # Check if resolution or number of images are within limits for radiation damage procedure
 ansR <- okRadDam(data)
 if (abs(fdrop) < 0.000001) ansR <- FALSE
 originalBatch <- max(data$V4,na.rm=TRUE)
 Ibatches <- c(Ibatches,range(data$V4,na.rm=TRUE)[1])
 Fbatches <- c(Fbatches,range(data$V4,na.rm=TRUE)[2])
 meanBC <- mean(data$V4,na.rm=TRUE)

 # Procedure for radiation damage detection only activated if resolution or number of images are within limits
 if (ansR)
 {
  wplots <- Wplots(data,nbin=nbin,m=0.125)
  tempWplots <- c(tempWplots,list(wplots))
  coefs <- newraddam(wplots)
 }
 if (!ansR) coefs <- c(1/0,1/0)

 # Purging is done only when criteria are met
 if (is.finite(coefs[2]))
 {
  cat(sprintf(" This datasets will be purged ..."))
  # Calculate image cutoff for each resolution
  tmp <- sapply(data$V9,tfit,coefs=coefs,fdrop=fdrop)+Ibatches[length(Ibatches)]
  tmp2 <- range(data$V4,na.rm=TRUE)
  tmp[tmp > tmp2[2]] <- NA
  tmp[tmp < 1] <- NA

  # Now purge data of radiation damaged ones
  cond <- data$V4 > tmp
  idx <- which(cond)
  if (length(idx) > 0) data <- data[-idx,]
  meanAC <- mean(data$V4,na.rm=TRUE)
  cat(sprintf(" done.\n"))
 }
 if (!is.finite(coefs[2]))
 {
  if (ansR) cat(sprintf(" This datasets will not be purged.\n"))
  if (!ansR) 
  {
   if (abs(fdrop) < 0.000001) cat(sprintf(" This datasets will not be purged.\n"))
   if (abs(fdrop) >= 0.000001) cat(sprintf(" This datasets will not be purged because of insufficient number of images or resolution too low.\n"))
  }
  meanAC <- meanBC
 }

 # Store purged datasets
 listData <- c(listData,list(data))

 # Store resolution ranges for later 
 sMax <- c(sMax,range(data$V9,na.rm=TRUE)[2])
 sMin <- c(sMin,range(data$V9,na.rm=TRUE)[1])

 # Number of last batch, after cut (approximate method at present)
 batches <- c(batches,floor(originalBatch-meanBC+meanAC))
}

#######################################################################################################
#######################################################################################################
###################################### DESCRIPTORS SECTION ############################################
######################################     B FACTORS       ############################################
#######################################################################################################
#######################################################################################################
# Find common resolution range for all datasets
m <- max(sMin)
M <- min(sMax)

# Compute overall wilson plots for all datasets
listBplots <- c()
listWplots <- c()
listSDWplots <- c()
for (data in listData)
{
 tmp <- aveIvsRes(nbin,data$V9^2,data$V7,m^2,M^2)
 listBplots <- c(listBplots,list(tmp))
 tmp <- aveIvsRes(nbin,data$V9,data$V7)
 listWplots <- c(listWplots,list(tmp))
 tmp <- aveIvsRes(nbin,data$V9,data$V8)
 listSDWplots <- c(listSDWplots,list(tmp))
}

# Get rid of listData to avoid dumping redundant information in R image
rm(listData)

# Choose reference dataset
tmp <- c()
for (wplots in listBplots) tmp <- c(tmp,sum(wplots[[1]],na.rm=TRUE))
idxref <- which(tmp == max(tmp))
idxref <- idxref[1]

#######################################################################################################
#######################################################################################################
########################################## RESOLUTION CUTS ############################################
#######################################################################################################
#######################################################################################################
# Signal-to-noise analysis to find out optimal resolution cut for each dataset
resoCuts <- c()
for (i in 1:length(listWplots))
{
 tmp <- bestReso(listWplots[[i]][[2]],listWplots[[i]][[3]]/listSDWplots[[i]][[3]],snr=isigi)
 resoCuts <- c(resoCuts,tmp[2])
 if (tmp[1]) 
 {
  nome <- filenames[i]
  tmp2 <- strsplit(nome,"_")
  tmp2 <- strsplit(tmp2[[1]][3],".",fixed=TRUE)
  idx <- as.integer(tmp2[[1]][1])
  messaggio <- sprintf("WARNING! Intensities of dataset %d are very weak. You might need to adjust the ISIGI parameter.\n",idx)
  cat(messaggio)
 }
}

# Get rid of listWplots and listSDWplots to avoid dumping redundant information in R image
rm(listWplots)
rm(listSDWplots)

#######################################################################################################
#######################################################################################################
###################################### DESCRIPTORS SECTION ############################################
######################################   CELL PARAMETERS   ############################################
#######################################################################################################
#######################################################################################################
# Load data from file "forR_macropar.dat"
macropar <- read.table("forR_macropar.dat")
names(macropar) <- c("cn","a","b","c","alpha","beta","gamma","mosa","ctoddist","wlength","centering")

# Re-arrange macropar with same order as file NEW_list_of_files.dat
macropar <- macropar[order(macropar$cn),]

# Prepare data frame for cluster analysis
cat("Preparing data for cluster analysis ...\n")
maindf <- cbind(macropar[,c(1,2,3,4,5,6,7,11)],data.frame(batch=batches,Ibatch=Ibatches,Fbatch=Fbatches))

# Decide which macroparameters can be used in cluster analysis
cell_columns <- c()
if (sum(is.na(maindf$a)) == 0)
{
 tmp <- sd(maindf$a)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,2)
}
if (sum(is.na(maindf$b)) == 0)
{
 tmp <- sd(maindf$b)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,3)
}
if (sum(is.na(maindf$c)) == 0)
{
 tmp <- sd(maindf$c)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,4)
}
if (sum(is.na(maindf$alpha)) == 0)
{
 tmp <- sd(maindf$alpha)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,5)
}
if (sum(is.na(maindf$beta)) == 0)
{
 tmp <- sd(maindf$beta)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,6)
}
if (sum(is.na(maindf$gamma)) == 0)
{
 tmp <- sd(maindf$gamma)
 if (tmp > 0.000001) cell_columns <- c(cell_columns,7)
}

#For NCdist, kill nParC
cell_columns <- c(1)


# Full list of descriptors used
fullc <- cell_columns

# Name lowresos and highresos arrays with same number as original dataset
lowresos <-  1/sMin
highresos <- resoCuts
names(lowresos) <- maindf$cn
names(highresos) <- maindf$cn

# New set of descriptors that better exploit information in Wilson plots
nparW <- nparWilson(listBplots,maindf$cn,nref=idxref,cumprob=0.95)

# Get rid of listBplots to avoid dumping redundant information in R image
rm(listBplots)

#######################################################################################################
#######################################################################################################
########################################## CLUSTER ANALYSIS ###########################################
#######################################################################################################
#######################################################################################################
# Enough information to start cluster analysis (add bit concerning cluster analysis here)
cat("Cluster analysis initiated ...\n")

# Carry on cluster analysis with descriptors from Wilson plots
nparW.dist <- dist(nparW)

# Make a dist object from the cell data
distMat<-evaluateSFChange(maindf)
rownames(distMat) <- maindf$cn
colnames(distMat) <- maindf$cn
distCAll <- as.dist(distMat)

# Final distance matrix has contributions from both set of descriptors
if (cparweight <= 0) distAll <- nparW.dist
if (cparweight > 0) distAll <- cparweight*distCAll+(1-cparweight)*nparW.dist

# Cluster analysis
tmp <- R.Version()
if (as.numeric(tmp[[6]]) >= 3 & as.numeric(tmp[[7]]) >= 1) {
 npar.hc_ward <- hclust(distAll,method="ward.D")
} else {
 npar.hc_ward <- hclust(distAll,method="ward")
}

# Dendrogram and resolution information are stored in list groups
groups <- treeToList(npar.hc_ward,lowresos,highresos)

# Calculate LCV values for all nodes of dendrogram
LCV_values <- c()
aLCV_values <- c()
LCV_couples <- matrix(nrow=length(groups[[1]]),ncol=2)
icpls <- 0
for (cn in groups[[1]])
{
 icpls <- icpls+1
 idx <- match(cn,macropar$cn)
 #tmp <- maxRatio(macropar[idx,2:7])
 tmp <- maxRatio(macropar,idx)
 LCV_values <- c(LCV_values,tmp[1])
 aLCV_values <- c(aLCV_values,tmp[2])
 LCV_couples[icpls,1] <- as.integer(tmp[3])
 LCV_couples[icpls,2] <- as.integer(tmp[4])
}

# Output merging nodes table to an ascii file
nTable <- "./CLUSTERS.txt"
if (file.exists(nTable)) emptyc <- file.remove(nTable)
linea <- "                               \n"
cat(linea,file=nTable)
linea <- " Cluster     Number of         Cluster         LCV      aLCV   Furthest    Datasets\n"
cat(linea,file=nTable,append=TRUE)
linea <- "  Number      Datasets          Height                         Datasets          ID\n"
cat(linea,file=nTable,append=TRUE)
linea <- "                               \n"
cat(linea,file=nTable,append=TRUE)
for (i in 1:length(groups[[1]]))
{
 sorted_groups <- sort(groups[[1]][[i]])
 linea <- paste(sprintf("     %03d           %3d         %7.3f     %7.2f %9.2f %5d %4d  ",
                i,length(groups[[1]][[i]]),npar.hc_ward$height[i],LCV_values[i],aLCV_values[i],LCV_couples[i,1],LCV_couples[i,2]),"  ",
                paste(sorted_groups,collapse=" "),"\n",sep="")
 #               i,length(groups[[1]][[i]]),npar.hc_ward$height[i],LCV_values[i],aLCV_values[i]),"  ",paste(groups[[1]][[i]],collapse=" "),"\n",sep="")
 cat(linea,file=nTable,append=TRUE)
}

# Print dendrogram with LCV values annotated (in PNG and POSTSCRIPT formats)
VV <- c()
for (i in 1:length(maindf[,1]))
{
 tmp <- .dcl(maindf[i,2],maindf[i,3],maindf[i,4],maindf[i,5],maindf[i,6],maindf[i,7])[16]
 attributes(tmp) <- NULL
 VV <- c(VV,tmp)
}
kk <- LCV_values[length(groups[[1]])]
akk <- aLCV_values[length(groups[[1]])]
tmsg <- sprintf("LCV: %.2f%s (absolute LCV: %.2f",kk,"%",akk)
msg <- bquote(.(tmsg) ~ ring(A)*")")

if (length(npar.hc_ward$height) > 1)
{
 png(file="./tree.png",height=1000,width=1000)
 #plclust(npar.hc_ward,xlab="Individual datasets",ylab="Ward distance",main=msg,sub="")
 par(oma=c(0,0.5,0,0.5),mgp=c(1.5,0.5,0))
 plot(npar.hc_ward,xlab="Individual datasets",ylab="Ward distance",main=msg,sub="",col.main="red",cex.main=2,col.lab="blue",cex.lab=2)
 nodesxy <- find_nodes_coords(npar.hc_ward,groups[[1]],macropar$cn)
 if (length(LCV_values) > 5) idx <- (length(LCV_values)-4):(length(LCV_values)-1)
 if (length(LCV_values) <= 5) idx <- 1:(length(LCV_values)-1)
 labelsxy <- c()
 for (i in 1:length(LCV_values))
 {
  stmp <- sprintf("%6.2f(%.2f)",LCV_values[i],aLCV_values[i])
  labelsxy <- c(labelsxy,stmp)
 }
 text(nodesxy$x[idx],nodesxy$y[idx],labels=labelsxy[idx],adj=c(0,-0.5),col="red",cex=1.0)
 emptyc <- dev.off()    # The "emptyc" is to collect the return value of dev.off(), so that it's not output
 #postscript(file="./tree.ps", height = 10, width = 10, paper = "a4")
 postscript(file="./tree.ps",paper = "a4",horizontal=FALSE)
 #plclust(npar.hc_ward,xlab="Individual datasets",ylab="Ward distance",main=msg,sub="")
 par(oma=c(0,0.5,0,0.5),mgp=c(1.5,0.5,0))
 plot(npar.hc_ward,xlab="Individual datasets",ylab="Ward distance",main=msg,sub="",col.main="red",cex.main=2,col.lab="blue",cex.lab=2)
 text(nodesxy$x[idx],nodesxy$y[idx],labels=labelsxy[idx],adj=c(0,-0.2),col="red",cex=0.8)
 emptyc <- dev.off()    # The "emptyc" is to collect the return value of dev.off(), so that it's not output
}
if (length(npar.hc_ward$height) == 1) cat("WARNING! No plot of dendrogram available as there is only 1 node. For cluster height refer to file 'CLUSTERS.txt'\n")
cat("Cluster analysis completed!\n")

# Remove "refs_*_*.dat" files and "forR_raddam.dat" file
for (crystal in filenames)
{
 if (file.exists(crystal)) emptyc <- file.remove(crystal)  # The "emptyc" is to collect the return value of file.remove(), so that it's not output
}
if (file.exists("forR_raddam.dat")) emptyc <- file.remove("forR_raddam.dat")
if (file.exists("forR_macropar.dat")) emptyc <- file.remove("forR_macropar.dat")

# Write file with batch information and original crystal number
filenames <- read.table("./NEW_list_of_files.dat",as.is=c(1))
emptyc <- file.remove("./NEW_list_of_files.dat")
for (ii in 1:length(filenames[,1]))
{
 stmp <- normalizePath(filenames$V1[ii])
 filenames$V1[ii] <- stmp
}
if (file.exists("FINAL_list_of_files.dat")) emptyc <- file.remove("FINAL_list_of_files.dat")
idx <- match(1:length(filenames[,1]),maindf$cn) 
for (i in 1:length(filenames[,1]))
{
 if (!is.na(idx[i]))
 {
  linea <- sprintf("%-130s %5.0f %5.0f %5.0f %5.0f %6.3f\n",filenames[i,],maindf$cn[idx[i]],
                   maindf$batch[idx[i]],maindf$Ibatch[idx[i]],maindf$Fbatch[idx[i]],highresos[idx[i]])
 }
 if (is.na(idx[i]))
 {
  linea <- sprintf("%-130s %5.0f %5.0f %5.0f %5.0f %6.3f\n",filenames[i,],NA,NA,NA,NA,NA)
 }
 cat(linea,file="./FINAL_list_of_files.dat",append=T)
}

# Output message for user
#cat("                     \n")
#cat("##################################################################\n")
#cat("##################################################################\n")
#cat("##################################################################\n")
#cat("   The following files have been produced:\n")
#cat("mtz_names.dat                 : ascii file containing path of input files used cwfor the analysis;\n")
#cat("BLEND_SUMMARY.txt             : ascii file containing a summary table with features of the datasets analysed;\n")
#cat("FINAL_list_of_files.dat       : ascii file with the path, serial number and number of batch used for each dataset;\n")
#cat("tree.png                      : PNG format picture of clustering dendrogram;\n")
#cat("CLUSTERS.txt                  : ascii file providing details of clustering.\n")
#cat("\n =========>>> Please, analyse the dendrogram and run BLEND in synthesis mode to produce merged datasets.\n")
#cat("    \n")

# Save image for next run of BLEND
save.image(file="BLEND.RData")

# Exit without saving
q(save = "no", status = 0, runLast = FALSE)
