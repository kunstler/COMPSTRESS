
####################################################################
####################################################################
####################################################################
###### CELLULAR SIMULATION MODEL 


## source functions
source(file="./R/RTheoModelFun.R") 




## plot fof how competition hierarchy model translate in probability of win between multiple species
par(mfrow=c(1,2))
plot((0:100)/100,function_compet(vec.t=(0:100)/100,K=1000),type="l",xlab="trait c",ylab="Probability of win") # ,log="y"
lines((0:100)/100,function_compet(vec.t=(0:100)/100,K=5),lty=3)
lines((0:100)/100,function_compet(vec.t=(0:100)/100,K=1),lty=4)
lines((0:100)/100,function_compet(vec.t=(0:100)/100,K=0.001),lty=6)

plot(seq(-1,1,length=100),1/(1+exp(-1000*(seq(-1,1,length=100)))),type="l",xlab="Diff in competitive ability c",ylab="Probability of wins")
lines(seq(-1,1,length=100),1/(1+exp(-5*(seq(-1,1,length=100)))),lty=3)
lines(seq(-1,1,length=100),1/(1+exp(-1*(seq(-1,1,length=100)))),lty=4)
lines(seq(-1,1,length=100),1/(1+exp(-0.01*(seq(-1,1,length=100)))),lty=6)



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
res.list.k100.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=100,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 5 P 0
res.list.k5.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 1 P 0
res.list.k1.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 0.001 P 0
res.list.k.001.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=0.001,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){


##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 100 P 1
res.list.k100.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=100,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 5 P 1
res.list.k5.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 1 P 1
res.list.k1.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 0.001 P 1
res.list.k.001.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=0.001,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

#########
## function run simulation
fun.run.sim <- function(A,B,Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

res.list.temp <- list()
res.list.temp[[1]] <- Alandscape.init
Alandscape <-  Alandscape.init

for(j in 1:A){
  for (i in 1:B){
  fun.update.landscape(fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=0.001,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)
  }  
res.list.temp[[j+1]] <- Alandscape
gc()
}
return(res.list.temp)
}




##### parameters to explore in a first set of simulation
## 
## fecundity  param.DISP=2
## compet hierarchy asymetry ,param.K=2 test k 100 5 1 0.001
##  premption  ,param.P=0

## need to think to the shape of the mortality function alternative option



length(res.list)

saveRDS(res.list.k100.p0,file="./output/res.list.k100.p0.rds")
saveRDS(res.list.k5.p0,file="./output/res.list.k5.p0.rds")
saveRDS(res.list.k1.p0,file="./output/res.list.k1.p0.rds")
saveRDS(res.list.k.001.p0,file="./output/res.list.k.001.p0.rds")
saveRDS(res.list.k100.p1,file="./output/res.list.k100.p1.rds")
saveRDS(res.list.k5.p1,file="./output/res.list.k5.p1.rds")
saveRDS(res.list.k1.p1,file="./output/res.list.k1.p1.rds")
saveRDS(res.list.k.001.p1,file="./output/res.list.k.001.p1.rds")

####read old RDS
res.list <- readRDS("./output/res.list.rds")


 ## plot landscape at different time step
par(mfrow=c(2,2))
image(x=1:nrow(res.list[[1]]),y=1:ncol(res.list[[1]]),z=res.list[[1]],xlim=c(1,ncol(res.list[[1]])),ylim=c(1,ncol(res.list[[1]])))
 image(x=1:nrow(res.list[[1]]),y=1:ncol(res.list[[1]]),z=res.list[[2]],xlim=c(1,ncol(res.list[[1]])),ylim=c(1,ncol(res.list[[1]])))
 image(x=1:nrow(res.list[[1]]),y=1:ncol(res.list[[1]]),z=res.list[[6]],xlim=c(1,ncol(res.list[[1]])),ylim=c(1,ncol(res.list[[1]])))
image(x=1:nrow(res.list[[1]]),y=1:ncol(res.list[[1]]),z=res.list[[11]],xlim=c(1,ncol(res.list[[1]])),ylim=c(1,ncol(res.list[[1]])))


par(mfcol=c(4,2),mar=c(2,2,2,2))
fun.plot.grad.quant(res.list.k100.p0)
fun.plot.grad.quant(res.list.k5.p0)
fun.plot.grad.quant(res.list.k1.p0)
fun.plot.grad.quant(res.list.k.001.p0)
fun.plot.grad.quant(res.list.k100.p1)
fun.plot.grad.quant(res.list.k5.p1)
fun.plot.grad.quant(res.list.k1.p1)
fun.plot.grad.quant(res.list.k.001.p1)


#### plot abundance
gray.col.vec <- rev(gray.colors(n=length(res.list)))
abun.temp <- fun.gradient.abundance.levels(res.list,t=1,imax=300)
plot(abun.temp[1,],type="l",ylim=c(0,1000),col=gray.col.vec[1])

for (i in 2:length(res.list)){
abun.temp <- fun.gradient.abundance.levels(res.list,t=i,imax=300)
lines(abun.temp[1,],col=gray.col.vec[i])
}

## proceesing output
