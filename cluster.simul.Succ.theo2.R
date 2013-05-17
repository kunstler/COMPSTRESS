
####################################################################
####################################################################
####################################################################
###### CELLULAR SIMULATION MODEL 

## K 100 5 1 0.001
# p 0
# disp 20

## source functions
source(file="./R/RTheoModelFun.R") 



##################################
## initialization  Succ
## Landscapoe matrix
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

##################################
##################################
#### RUN A SIMULATION SUCCESSION TES  FOR A B time step
## K 100 P 1 Succ 0.5
res.list.Succ.k100.p0.Succ0.5 <- fun.run.sim.Succ(A=20,B=50,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=100,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.5,dist.vec,array.i,array.j)

saveRDS(res.list.Succ.k100.p0.Succ0.5,file="./output/res.list.Succ.k100.p0.Succ0.5.rds")

##################################
##################################
#### RUN A SIMULATION SUCCESSION TES  FOR A B time step
## K 5 P 1 Succ 0.5
res.list.Succ.k5.p0.Succ0.5 <- fun.run.sim.Succ(A=20,B=50,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.5,dist.vec,array.i,array.j)

saveRDS(res.list.Succ.k5.p0.Succ0.5,file="./output/res.list.Succ.k5.p0.Succ0.5.rds")

##################################
##################################
#### RUN A SIMULATION SUCCESSION TES  FOR A B time step
## K 1 P 1 Succ 0.5
res.list.Succ.k1.p0.Succ0.5 <- fun.run.sim.Succ(A=20,B=50,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.5,dist.vec,array.i,array.j)

saveRDS(res.list.Succ.k1.p0.Succ0.5,file="./output/res.list.Succ.k1.p0.Succ0.5.rds")

##################################
##################################
#### RUN A SIMULATION SUCCESSION TES  FOR A B time step
## K 0.001 P 1 Succ 0.5
res.list.Succ.k.001.p0.Succ0.5 <- fun.run.sim.Succ(A=20,B=50,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=0.001,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.5,dist.vec,array.i,array.j)


saveRDS(res.list.Succ.k.001.p0.Succ0.5,file="./output/res.list.Succ.k.001.p0.Succ0.5.rds")
