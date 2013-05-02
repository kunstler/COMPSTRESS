
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
## initialization
## Landscapoe matrix
Nlandscape <-  3 #dimension (in 2*NN) of the long side corressponding to the climatic stress gradient
NN <- 100  ## size of landscape (classicaly 256)
N <-  1 ## number of neigborhood cell 1 = moor neigborhood
Alandscape <-  matrix(NA,nrow=NN,ncol=Nlandscape*NN)
## climate gradient
climate.grad <-  seq(from=0,to=1,length=Nlandscape*NN)

## init landscape with random draw of competition traits
init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*0.5))
vec.land <- as.vector(Alandscape)
vec.land[init.temp] <- sample((0:100)/100,size=round(NN*NN*Nlandscape*0.5),replace=T)
Alandscape <- matrix(vec.land,nrow=NN,ncol=NN*Nlandscape)
Alandscape.init <-  Alandscape

## create a table of 8 neigborhood cells

list.temp <- function.array.neigcells(NN,Nlandscape)

array.i <- list.temp[[1]]
array.j <- list.temp[[2]]
rm(list.temp)
array.i[10,20,]



##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 100 P 0
res.list.k100.p0.disp20 <- fun.run.sim(A=20,B=50,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=20,param.K=100,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 5 P 0
res.list.k5.p0.disp20 <- fun.run.sim(A=20,B=50,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=20,param.K=5,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 1 P 0
res.list.k1.p0.disp20 <- fun.run.sim(A=20,B=50,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=20,param.K=1,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 0.001 P 0
res.list.k.001.p0.disp20 <- fun.run.sim(A=20,B=50,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=20,param.K=0.001,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)


saveRDS(res.list.k100.p0.disp20,file="./output/res.list.k100.p0.disp20.rds")
saveRDS(res.list.k5.p0.disp20,file="./output/res.list.k5.p0.disp20.rds")
saveRDS(res.list.k1.p0.disp20,file="./output/res.list.k1.p0.disp20.rds")
saveRDS(res.list.k.001.p0.disp20,file="./output/res.list.k.001.p0.disp20.rds")
