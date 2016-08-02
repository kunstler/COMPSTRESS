####################################################################
####################################################################
####################################################################
###### Test CELLULAR SIMULATION MODEL in Rcpp to speed up

## source functions
library(Rcpp)
sourceCpp(file.path('R','Rcpp', 'CS_Cell.cpp'))


## Initialize the data in R (to be move in a function)
Nlandscape <-  3 #dimension (in 2*NN) of the long side corressponding to the climatic stress gradient
NN <- 100  ## size of landscape (classicaly 256)
N <-  1 ## number of neigborhood cell 1 = moor neigborhood
AlandscapeE <-  matrix(NA,nrow=NN,ncol=Nlandscape*NN)
AlandscapeL <- AlandscapeE
Alandscape.Succ <- matrix("NO",nrow=NN,ncol=Nlandscape*NN)

## climate gradient
climate.grad <-  seq(from=0,to=1,length=Nlandscape*NN)

## init landscape with random draw of competition traits
init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*0.5))
vec.land <- as.vector(AlandscapeE)
vec.landE <- as.vector(AlandscapeE)
vec.landL <- as.vector(AlandscapeE)
vec.land.Succ <- as.vector(Alandscape.Succ)

### create random species with a tradeoff between competition stress tolerance and eraly and late succ compet
vec.land[init.temp] <- sample((0:100),size=round(NN*NN*Nlandscape*0.5),replace=T) ## stress tolerance

## fun to create random early succ and late suss competitive ability
fun.sample.compet <-  function(x){
  sample(0:x,size=1)
}
vec.landE[init.temp] <- sapply(vec.land[init.temp],FUN=fun.sample.compet)
vec.landL[init.temp] <-(vec.land[init.temp] -vec.landE[init.temp])/100
vec.landE[init.temp] <- vec.landE[init.temp]/100

## par(mfrow=c(2,2))
## x   <- cbind(vec.landE[init.temp] ,vec.landL[init.temp])
## data <- as.data.frame(x)
## names(data) <- c("X","Y")
## colors  <- densCols(x)
## plot(x, col=colors, pch=20,cex=0.25,xlab="C.E",ylab="C.L")
## x   <- cbind(1-vec.land[init.temp] ,vec.landL[init.temp])
## data <- as.data.frame(x)
## names(data) <- c("X","Y")
## colors  <- densCols(x)
## plot(x, col=colors, pch=20,cex=0.25,xlab="S",ylab="C.L")
## x   <- cbind(vec.landE[init.temp] ,1-vec.land[init.temp])
## data <- as.data.frame(x)
## names(data) <- c("X","Y")
## colors  <- densCols(x)
## plot(x, col=colors, pch=20,cex=0.25,xlab="C.E",ylab="S")


vec.land.Succ[init.temp] <-  rep("E",round(NN*NN*Nlandscape*0.5))
AlandscapeE <- matrix(vec.landE,nrow=NN,ncol=NN*Nlandscape)
AlandscapeL <- matrix(vec.landL,nrow=NN,ncol=NN*Nlandscape)
Alandscape.Succ <- matrix(vec.land.Succ,nrow=NN,ncol=NN*Nlandscape)

Alandscape.LIST.init <-  list(AlandscapeE,AlandscapeL,Alandscape.Succ)


## create a table of 8 neigborhood cells

list.temp <- function.array.neigcells(NN,Nlandscape)

array.i <- list.temp[[1]]
array.j <- list.temp[[2]]
rm(list.temp)

