
####################################################################
####################################################################
####################################################################
###### CELLULAR SIMULATION MODEL 


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


######
## OTHER PARAMETERS TO EXPLORE

## param.DISP=2
### param.dist=0.1
## rather explore the ratio of the two 2/0.1

## test param.DISP=20 or 0.2

### other trade off ??
# shape of moratlity function
# add constant competition from inferior
# multiple trade off
# colonization ? 

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 100 P 0
res.list.k100.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=100,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 5 P 0
res.list.k5.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 1 P 0
res.list.k1.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 0.001 P 0
res.list.k.001.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=0.001,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)


##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 100 P 1
res.list.k100.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=100,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 5 P 1
res.list.k5.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 1 P 1
res.list.k1.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 0.001 P 1
res.list.k.001.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=0.001,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)




##### parameters to explore in a first set of simulation
## 
## fecundity  param.DISP=2
## compet hierarchy asymetry ,param.K=2 test k 100 5 1 0.001
##  premption  ,param.P=0

## need to think to the shape of the mortality function alternative option

saveRDS(res.list.k100.p0,file="./output/res.list.k100.p0.rds")
saveRDS(res.list.k5.p0,file="./output/res.list.k5.p0.rds")
saveRDS(res.list.k1.p0,file="./output/res.list.k1.p0.rds")
saveRDS(res.list.k.001.p0,file="./output/res.list.k.001.p0.rds")
saveRDS(res.list.k100.p1,file="./output/res.list.k100.p1.rds")
saveRDS(res.list.k5.p1,file="./output/res.list.k5.p1.rds")
saveRDS(res.list.k1.p1,file="./output/res.list.k1.p1.rds")
saveRDS(res.list.k.001.p1,file="./output/res.list.k.001.p1.rds")

####read old RDS
res.list.k100.p0 <- readRDS(file="./output/res.list.k100.p0.rds")
res.list.k5.p0 <- readRDS(file="./output/res.list.k5.p0.rds")
res.list.k1.p0 <- readRDS(file="./output/res.list.k1.p0.rds")
res.list.k.001.p0 <- readRDS(file="./output/res.list.k.001.p0.rds")
res.list.k100.p1 <- readRDS(file="./output/res.list.k100.p1.rds")
res.list.k5.p1 <- readRDS(file="./output/res.list.k5.p1.rds")
res.list.k1.p1 <- readRDS(file="./output/res.list.k1.p1.rds")
res.list.k.001.p1 <- readRDS(file="./output/res.list.k.001.p1.rds")

 ## plot landscape at different time step ,mar=c(3,2,2,3,2)
par(mfrow=c(2,4),mar=c(1,4,1,1),oma=c(1,4,2,1))
image(x=1:nrow(res.list.k100.p0[[15]]),y=1:ncol(res.list.k100.p0[[15]]),z=res.list.k100.p0[[15]],asp=1,xlab=NA,ylab="cliamte gradient")
image(x=1:nrow(res.list.k5.p0[[15]]),y=1:ncol(res.list.k5.p0[[15]]),z=res.list.k5.p0[[15]],asp=1,xlab=NA,ylab="cliamte gradient")
image(x=1:nrow(res.list.k1.p0[[15]]),y=1:ncol(res.list.k1.p0[[15]]),z=res.list.k1.p0[[15]],asp=1,xlab=NA,ylab="cliamte gradient")
image(x=1:nrow(res.list.k.001.p0[[15]]),y=1:ncol(res.list.k.001.p0[[15]]),z=res.list.k.001.p0[[15]],asp=1,xlab=NA,ylab="cliamte gradient")
#
image(x=1:nrow(res.list.k100.p1[[15]]),y=1:ncol(res.list.k100.p1[[15]]),z=res.list.k100.p1[[15]],asp=1,xlab=NA,ylab="cliamte gradient")
image(x=1:nrow(res.list.k5.p1[[15]]),y=1:ncol(res.list.k5.p1[[15]]),z=res.list.k5.p1[[15]],asp=1,xlab=NA,ylab="cliamte gradient")
image(x=1:nrow(res.list.k1.p1[[15]]),y=1:ncol(res.list.k1.p1[[15]]),z=res.list.k1.p1[[15]],asp=1,xlab=NA,ylab="cliamte gradient")
image(x=1:nrow(res.list.k.001.p1[[15]]),y=1:ncol(res.list.k.001.p1[[15]]),z=res.list.k.001.p1[[15]],asp=1,xlab=NA,ylab="cliamte gradient")

## labels
mtext("Premption", 2,  adj=0.80,padj=-1.5, cex=1.2, outer=TRUE)
mtext("No Premption",2,adj=0.2,padj=-1.5, cex=1.2, outer=TRUE)

 mtext("K=100", 3   , adj=0.1, cex=1.2, outer=TRUE)
 mtext("K=5",  3   , adj=0.37, cex=1.2, outer=TRUE)
 mtext("K=100", 3  , adj=0.65, cex=1.2, outer=TRUE)
 mtext("K=0.001", 3   , adj=0.92, cex=1.2, outer=TRUE)



par(mfcol=c(4,2),mar=c(4,4,2,2))
fun.plot.grad.quant(res.list.k100.p0)
fun.plot.grad.quant(res.list.k5.p0)
fun.plot.grad.quant(res.list.k1.p0)
fun.plot.grad.quant(res.list.k.001.p0)
#
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


