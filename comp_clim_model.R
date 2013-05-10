
####################################################################
####################################################################
####################################################################
###### CELLULAR SIMULATION MODEL 


## source functions
source(file="./R/RTheoModelFun.R") 



##################################
## initialization NO Succ
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
res.list.k100.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=100,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)



##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 5 P 0
res.list.k5.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 1 P 0
res.list.k1.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
 disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 0.001 P 0
res.list.k.001.p0 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
 disp.fun=disp.unif.fun,param.DISP=2,param.K=0.001,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)


##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 100 P 1
res.list.k100.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=100,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 5 P 1
res.list.k5.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
 disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 1 P 1
res.list.k1.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

##################################
##################################
#### RUN A SIMULATION FOR A B time step
## K 0.001 P 1
res.list.k.001.p1 <- fun.run.sim(A=15,B=20,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=0.001,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)




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


##################################3
##################################
### TEST SUCCESSION

##################################
## initialization  Succ
## Landscapoe matrix
Nlandscape <-  3 #dimension (in 2*NN) of the long side corressponding to the climatic stress gradient
NN <- 100  ## size of landscape (classicaly 256)
N <-  1 ## number of neigborhood cell 1 = moor neigborhood
AlandscapeE <-  matrix(NA,nrow=NN,ncol=Nlandscape*NN)
AlandscapeL <- AlandscapeE
Alandscape.Succ <- AlandscapeE

## climate gradient
climate.grad <-  seq(from=0,to=1,length=Nlandscape*NN)

## init landscape with random draw of competition traits
init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*0.5))
vec.land <- as.vector(Alandscape)
vec.landE <- as.vector(Alandscape)
vec.landL <- as.vector(Alandscape)
vec.land.Succ <- as.vector(Alandscape)
vec.land[init.temp] <- sample((0:100)/100,size=round(NN*NN*Nlandscape*0.5),replace=T)
vec.landE[init.temp] <- vec.land[init.temp]*sample((0:100)/100,size=round(NN*NN*Nlandscape*0.5),replace=T)
vec.landL[init.temp] <-vec.land[init.temp] -vec.landE[init.temp]
vec.land.Succ[init.temp] <-  rep("E",round(NN*NN*Nlandscape*0.5))
AlandscapeE <- matrix(vec.landE,nrow=NN,ncol=NN*Nlandscape)
AlandscapeL <- matrix(vec.landL,nrow=NN,ncol=NN*Nlandscape)
Alandscape.Succ <- matrix(vec.land.Succ,nrow=NN,ncol=NN*Nlandscape)
AlandscapeE.init <-  AlandscapeE
AlandscapeL.init <-  AlandscapeL
Alandscape.Succ.init <-  Alandscape.Succ

## create a table of 8 neigborhood cells

list.temp <- function.array.neigcells(NN,Nlandscape)

array.i <- list.temp[[1]]
array.j <- list.temp[[2]]
rm(list.temp)
array.i[10,20,]




##################################
##################################
#### RUN A SIMULATION SUCCESSION TES  FOR A B time step
## K 5 P 1
res.list.k.001.p1 <- fun.run.sim.Succ(A=2,B=2,Alandscape.init=Alandscape.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.25,dist.vec,array.i,array.j)




### NEED to read all outputs from cluster.

## read file in output cluster
list.files(path="./output_cluster",pattern="p0.disp20.r")

### create vector names of files to read from cluster for the plots
names.p0.disp20.dist0.1 <- c( "res.list.k100.p0.disp20.rds" ,  "res.list.k5.p0.disp20.rds", "res.list.k1.p0.disp20.rds" ,"res.list.k.001.p0.disp20.rds" )
names.p0.5.disp20.dist0.1 <- c( "res.list.k100.p0.5.disp20.rds" ,  "res.list.k5.p0.5.disp20.rds", "res.list.k1.p0.5.disp20.rds" ,"res.list.k.001.p0.5.disp20.rds" )
names.p1.disp20.dist0.1 <- c( "res.list.k100.p1.disp20.rds" ,  "res.list.k5.p1.disp20.rds", "res.list.k1.p1.disp20.rds" ,"res.list.k.001.p1.disp20.rds" )

names.p0.disp20.dist0.5 <- c( "res.list.k100.p0.disp20.dist0.5.rds" ,  "res.list.k5.p0.disp20.dist0.5.rds", "res.list.k1.p0.disp20.dist0.5.rds" ,"res.list.k.001.p0.disp20.dist0.5.rds" )
names.p0.5.disp20.dist0.5 <- c( "res.list.k100.p0.5.disp20.dist0.5.rds" ,  "res.list.k5.p0.5.disp20.dist0.5.rds", "res.list.k1.p0.5.disp20.dist0.5.rds" ,"res.list.k.001.p0.5.disp20.dist0.5.rds" )
names.p1.disp20.dist0.5 <- c( "res.list.k100.p1.disp20.dist0.5.rds" ,  "res.list.k5.p1.disp20.dist0.5.rds", "res.list.k1.p1.disp20.dist0.5.rds" ,"res.list.k.001.p1.disp20.dist0.5.rds" )


names.p0.disp2.dist0.1 <- c( "res.list.k100.p0.disp2.rds" ,  "res.list.k5.p0.disp2.rds", "res.list.k1.p0.disp2.rds" ,"res.list.k.001.p0.disp2.rds" )
names.p0.5.disp2.dist0.1 <- c( "res.list.k100.p0.5.disp2.rds" ,  "res.list.k5.p0.5.disp2.rds", "res.list.k1.p0.5.disp2.rds" ,"res.list.k.001.p0.5.disp2.rds" )
names.p1.disp2.dist0.1 <- c( "res.list.k100.p1.disp2.rds" ,  "res.list.k5.p1.disp2.rds", "res.list.k1.p1.disp2.rds" ,"res.list.k.001.p1.disp2.rds" )

names.p0.disp.2.dist0.1 <- c( "res.list.k100.p0.disp.2.rds" ,  "res.list.k5.p0.disp.2.rds", "res.list.k1.p0.disp.2.rds" ,"res.list.k.001.p0.disp.2.rds" )
names.p0.5.disp.2.dist0.1 <- c( "res.list.k100.p0.5.disp.2.rds" ,  "res.list.k5.p0.5.disp.2.rds", "res.list.k1.p0.5.disp.2.rds" ,"res.list.k.001.p0.5.disp.2.rds" )
names.p1.disp.2.dist0.1 <- c( "res.list.k100.p1.disp.2.rds" ,  "res.list.k5.p1.disp.2.rds", "res.list.k1.p1.disp.2.rds" ,"res.list.k.001.p1.disp.2.rds" )


#################
### PLOTS OUTPUTS OF CLUSTER SIMULATION


## disp20
pdf("./figs/Firstsimul.F20.pdf")
par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(3,4,2,1))
lapply(names.p0.disp20.dist0.1 ,fun.plot.grad.quant.cluster,
       path="output_cluster")
lapply(names.p0.5.disp20.dist0.1 ,fun.plot.grad.quant.cluster,
       path="output_cluster")
lapply(names.p1.disp20.dist0.1 ,fun.plot.grad.quant.cluster,
       path="output_cluster")
mtext("No Premption (P=0)",3,adj=0.1, cex=0.9, outer=TRUE)
mtext("Half  Premption (P=0.5)",3,adj=0.5, cex=0.9, outer=TRUE)
mtext("Premption (P=1)",3,adj=0.9, cex=0.9, outer=TRUE)
mtext("K=100", 2   , adj=0.92,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=5",  2   , adj=0.65,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=1", 2  , adj=0.37,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=0.001", 2   , adj=0.1,padj=-1.5, cex=1.2, outer=TRUE)
mtext("Climatic gradient", 1   , adj=0.5,padj=1.3, cex=1.4, outer=TRUE)
dev.off()

## par(mfcol=c(4,2),mar=c(4,4,2,2))
## lapply(names.p0.disp20.dist0.5 ,fun.plot.grad.quant.cluster,path="output_cluster")
## lapply(names.p1.disp20.dist0.5 ,fun.plot.grad.quant.cluster,path="output_cluster")

## disp2
pdf("./figs/Firstsimul.F2.pdf")
par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(3,4,2,1))
lapply(names.p0.disp2.dist0.1 ,fun.plot.grad.quant.cluster
       ,path="output_cluster")
lapply(names.p0.5.disp2.dist0.1 ,fun.plot.grad.quant.cluster
       ,path="output_cluster")
lapply(names.p1.disp2.dist0.1 ,fun.plot.grad.quant.cluster
       ,path="output_cluster")
mtext("No Premption (P=0)",3,adj=0.1, cex=0.9, outer=TRUE)
mtext("Half  Premption (P=0.5)",3,adj=0.5, cex=0.9, outer=TRUE)
mtext("Premption (P=1)",3,adj=0.9, cex=0.9, outer=TRUE)
mtext("K=100", 2   , adj=0.92,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=5",  2   , adj=0.65,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=1", 2  , adj=0.37,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=0.001", 2   , adj=0.1,padj=-1.5, cex=1.2, outer=TRUE)
mtext("Climatic gradient", 1   , adj=0.5,padj=1.3, cex=1.4, outer=TRUE)
dev.off()

## disp 0.2
pdf("./figs/Firstsimul.F.2.pdf")
par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(3,4,2,1))
lapply(names.p0.disp.2.dist0.1 ,fun.plot.grad.quant.cluster,
       path="output_cluster")
lapply(names.p0.5.disp.2.dist0.1 ,fun.plot.grad.quant.cluster,
       path="output_cluster")
lapply(names.p1.disp.2.dist0.1 ,fun.plot.grad.quant.cluster,
       path="output_cluster")
mtext("No Premption (P=0)",3,adj=0.1, cex=0.9, outer=TRUE)
mtext("Half  Premption (P=0.5)",3,adj=0.5, cex=0.9, outer=TRUE)
mtext("Premption (P=1)",3,adj=0.9, cex=0.9, outer=TRUE)
mtext("K=100", 2   , adj=0.92,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=5",  2   , adj=0.65,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=1", 2  , adj=0.37,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=0.001", 2   , adj=0.1,padj=-1.5, cex=1.2, outer=TRUE)
mtext("Climatic gradient", 1   , adj=0.5,padj=1.3, cex=1.4, outer=TRUE)
dev.off()


################### image

## plot landscape at last time step DISP20
pdf("./figs/Firstsimul.image.F20.pdf")
par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(3,4,2,1))
lapply(names.p0.disp20.dist0.1 ,fun.image.landscape,
       path="output_cluster")
lapply(names.p0.5.disp20.dist0.1 ,fun.image.landscape,
       path="output_cluster")
lapply(names.p1.disp20.dist0.1 ,fun.image.landscape,
       path="output_cluster")
mtext("No Premption (P=0)",3,adj=0.1, cex=0.9, outer=TRUE)
mtext("Half  Premption (P=0.5)",3,adj=0.5, cex=0.9, outer=TRUE)
mtext("Premption (P=1)",3,adj=0.9, cex=0.9, outer=TRUE)
mtext("K=100", 2   , adj=0.92,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=5",  2   , adj=0.65,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=1", 2  , adj=0.37,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=0.001", 2   , adj=0.1,padj=-1.5, cex=1.2, outer=TRUE)
dev.off()

## plot landscape at last time step DISP2
pdf("./figs/Firstsimul.image.F2.pdf")
par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(3,4,2,1))
lapply(names.p0.disp2.dist0.1 ,fun.image.landscape,
       path="output_cluster")
lapply(names.p0.5.disp2.dist0.1 ,fun.image.landscape,
       path="output_cluster")
lapply(names.p1.disp2.dist0.1 ,fun.image.landscape,
       path="output_cluster")
mtext("No Premption (P=0)",3,adj=0.1, cex=0.9, outer=TRUE)
mtext("Half  Premption (P=0.5)",3,adj=0.5, cex=0.9, outer=TRUE)
mtext("Premption (P=1)",3,adj=0.9, cex=0.9, outer=TRUE)
mtext("K=100", 2   , adj=0.92,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=5",  2   , adj=0.65,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=1", 2  , adj=0.37,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=0.001", 2   , adj=0.1,padj=-1.5, cex=1.2, outer=TRUE)
dev.off()

#### plot abundance
gray.col.vec <- rev(gray.colors(n=length(res.list)))
abun.temp <- fun.gradient.abundance.levels(res.list,t=1,imax=300)
plot(abun.temp[1,],type="l",ylim=c(0,1000),col=gray.col.vec[1])

for (i in 2:length(res.list)){
abun.temp <- fun.gradient.abundance.levels(res.list,t=i,imax=300)
lines(abun.temp[1,],col=gray.col.vec[i])
}


