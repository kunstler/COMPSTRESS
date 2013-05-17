
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
NN <- 20  ## size of landscape (classicaly 256)
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

## plot 3d trade off 
fun.trade.off.3d <- function(X,Y){
   c <- (1-(X+Y))

   if(sum(c<0)) {c[c<0] <- NA}
return(c)}

fun.trade.off.3d.1 <- function(X,Y){
 c0 <- (1-(X+Y))

  c <- ((X))
   if(sum(c0<0)) {c[c0<0] <- NA}
   return(c)}

fun.trade.off.3d.2 <- function(X,Y){
 c0 <- (1-(X+Y))
   c <- (Y)
   if(sum(c0<0)) {c[c0<0] <- NA}
   return(c)}

fun.trade.off.3d.3 <- function(X,Y){
   c <- (1-(X+Y))
   if(sum(c<0)) {c[c<0] <- NA}
   return(c)}

### TWO PLOTS
s.matrix <- outer(X=(0:100)/100,Y=(0:100)/100,FUN=fun.trade.off.3d) 
pdf("./figs/3dtrade.off.pdf")
par(mfrow=c(1,2))
image(x=(0:100)/100,y=(0:100)/100,z=s.matrix,xlab="Early Successional Competitive Ability",ylab="Late Successional Competitive Ability",col = terrain.colors(100))
dev.off()

 trade.array <- stack(flip(raster(outer(X=(0:100)/100,Y=(0:100)/100,FUN=fun.trade.off.3d.2) ),2),
                     flip(raster(outer(X=(0:100)/100,Y=(0:100)/100,FUN=fun.trade.off.3d.1)),2),
                     flip(raster(outer(X=(0:100)/100,Y=(0:100)/100,FUN=fun.trade.off.3d.3)),2))

pdf("./figs/3dtrade.offRGB.pdf")

plotRGB(trade.array,scale=1,axes=T,xlab="Early Successional Competitive Ability",ylab="Late Successional Competitive Ability")
text(0.15,0.1,labels="Climate stress tolerant")
text(0.9,0.1,labels="Early successional")
text(0.15,0.9,labels="Late successional")
dev.off()


##################################
## initialization  Succ
## Landscapoe matrix
Nlandscape <-  3 #dimension (in 2*NN) of the long side corressponding to the climatic stress gradient
NN <- 20  ## size of landscape (classicaly 256)
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

plot(vec.landE[init.temp] ,vec.landE[init.temp] + vec.landL[init.temp]  )
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
array.i[10,20,]




##################################
##################################
#### RUN A SIMULATION SUCCESSION TEST  FOR A B time step
## K 5 P 1
res.list.Succ.k5.p1.Succ0.95 <- fun.run.sim.Succ(A=20,B=5,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.15,dist.vec,array.i,array.j)

par(mfrow=c(5,5),mar=c(0,0,0,0))
for (i in c(0:20)) {image(res.list.Succ.k5.p1.Succ0.95[[1]][[i+1]] )
                    print(mean(res.list.Succ.k5.p1.Succ0.95[[1]][[i+1]],na.rm=T))}
x11()
par(mfrow=c(5,5),mar=c(0,0,0,0))
for (i in c(0:20)) {image(res.list.Succ.k5.p1.Succ0.95[[2]][[i+1]] )
                    print(mean(res.list.Succ.k5.p1.Succ0.95[[2]][[i+1]],na.rm=T))}


############ TEST IF SUCCESSION MODEL ALLOWS COEXISTENCE
##################################
## initialization  Succ
## Landscapoe matrix
Nlandscape <-  1 #dimension (in 2*NN) of the long side corressponding to the climatic stress gradient
NN <- 80  ## size of landscape (classicaly 256)
N <-  1 ## number of neigborhood cell 1 = moor neigborhood
AlandscapeE <-  matrix(NA,nrow=NN,ncol=Nlandscape*NN)
AlandscapeL <- AlandscapeE
Alandscape.Succ <- matrix("NO",nrow=NN,ncol=Nlandscape*NN)

## climate gradient
climate.grad <-  rep(0,length=Nlandscape*NN)

## init landscape with random draw of competition traits
init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*0.5))
vec.land <- as.vector(AlandscapeE)
vec.landE <- as.vector(AlandscapeE)
vec.landL <- as.vector(AlandscapeE)
vec.land.Succ <- as.vector(Alandscape.Succ)

### create random species with a tradeoff between competition stress tolerance and eraly and late succ compet
vec.land[init.temp] <- rep(100,length=round(NN*NN*Nlandscape*0.5)) ## stress tolerance

## fun to create random early succ and late suss competitive ability
fun.sample.compet <-  function(x){
  sample(0:x,size=1)
}
vec.landE[init.temp] <- sapply(vec.land[init.temp],FUN=fun.sample.compet)
vec.landL[init.temp] <-(vec.land[init.temp] -vec.landE[init.temp])/100
vec.landE[init.temp] <- vec.landE[init.temp]/100

plot(vec.landE[init.temp] ,vec.landL[init.temp]  )


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
array.i[10,20,]

res.list.Succ.k5.p1.Succ0.5 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.5,dist.vec,array.i,array.j)


pdf("./figs/Succ0.5.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.5[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.5[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.5[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.5[[2]][[20+1]])
dev.off()

# Succ 0.25
res.list.Succ.k5.p1.Succ0.25 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.25,dist.vec,array.i,array.j)


pdf("./figs/Succ0.25.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.25[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.25[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.25[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.25[[2]][[20+1]])
dev.off()

# Succ 0.35
res.list.Succ.k5.p1.Succ0.35 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.35,dist.vec,array.i,array.j)


pdf("./figs/Succ0.35.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.35[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.35[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.35[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.35[[2]][[20+1]])
dev.off()

# Succ 0.85
res.list.Succ.k5.p1.Succ0.85 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.85,dist.vec,array.i,array.j)


pdf("./figs/Succ0.85.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.85[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.85[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.85[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.85[[2]][[20+1]])
dev.off()

# Succ 0.05
res.list.Succ.k5.p1.Succ0.05 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.05,dist.vec,array.i,array.j)

pdf("./figs/Succ0.05.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.05[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.05[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.05[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.05[[2]][[20+1]])
dev.off()

# Succ 0.15
res.list.Succ.k5.p1.Succ0.15 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.15,dist.vec,array.i,array.j)

pdf("./figs/Succ0.15.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.15[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.15[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.15[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.15[[2]][[20+1]])
dev.off()

# Succ 0.20
res.list.Succ.k5.p1.Succ0.20 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.20,dist.vec,array.i,array.j)

pdf("./figs/Succ0.20.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.20[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.20[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.20[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.20[[2]][[20+1]])
dev.off()

# Succ 0.30
res.list.Succ.k5.p1.Succ0.30 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.30,dist.vec,array.i,array.j)

pdf("./figs/Succ0.30.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.30[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.30[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.30[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.30[[2]][[20+1]])
dev.off()

# Succ 0.95
res.list.Succ.k5.p1.Succ0.95 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.95,dist.vec,array.i,array.j)

pdf("./figs/Succ0.95.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.95[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.95[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.95[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.95[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.95,file="./output/res.list.Succ.k5.p1.Succ0.95.rds")
saveRDS(res.list.Succ.k5.p1.Succ0.85,file="./output/res.list.Succ.k5.p1.Succ0.85.rds")
saveRDS(res.list.Succ.k5.p1.Succ0.5,file="./output/res.list.Succ.k5.p1.Succ0.5.rds")
saveRDS(res.list.Succ.k5.p1.Succ0.35,file="./output/res.list.Succ.k5.p1.Succ0.35.rds")
saveRDS(res.list.Succ.k5.p1.Succ0.25,file="./output/res.list.Succ.k5.p1.Succ0.25.rds")
saveRDS(res.list.Succ.k5.p1.Succ0.15,file="./output/res.list.Succ.k5.p1.Succ0.15.rds")
saveRDS(res.list.Succ.k5.p1.Succ0.05,file="./output/res.list.Succ.k5.p1.Succ0.05.rds")

# Succ 0.24
res.list.Succ.k5.p1.Succ0.24 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.24,dist.vec,array.i,array.j)

pdf("./figs/Succ0.24.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.24[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.24[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.24[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.24[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.24,file="./output/res.list.Succ.k5.p1.Succ0.24.rds")
# Succ 0.23
res.list.Succ.k5.p1.Succ0.23 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.23,dist.vec,array.i,array.j)

pdf("./figs/Succ0.23.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.23[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.23[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.23[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.23[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.23,file="./output/res.list.Succ.k5.p1.Succ0.23.rds")

# Succ 0.22
res.list.Succ.k5.p1.Succ0.22 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.22,dist.vec,array.i,array.j)

pdf("./figs/Succ0.22.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.22[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.22[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.22[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.22[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.22,file="./output/res.list.Succ.k5.p1.Succ0.22.rds")

# Succ 0.21
res.list.Succ.k5.p1.Succ0.21 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.21,dist.vec,array.i,array.j)

pdf("./figs/Succ0.21.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.21[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.21[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.21[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.21[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.21,file="./output/res.list.Succ.k5.p1.Succ0.21.rds")

# Succ 0.20
res.list.Succ.k5.p1.Succ0.20 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.20,dist.vec,array.i,array.j)

pdf("./figs/Succ0.20.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.20[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.20[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.20[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.20[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.20,file="./output/res.list.Succ.k5.p1.Succ0.20.rds")

# Succ 0.1
res.list.Succ.k5.p1.Succ0.1 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.1,dist.vec,array.i,array.j)

pdf("./figs/Succ0.1.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.1[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.1[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.1[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.1[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.1,file="./output/res.list.Succ.k5.p1.Succ0.1.rds")

# Succ 0.30
res.list.Succ.k5.p1.Succ0.30 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.30,dist.vec,array.i,array.j)

pdf("./figs/Succ0.30.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.30[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.30[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.30[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.30[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.30,file="./output/res.list.Succ.k5.p1.Succ0.30.rds")

# Succ 0.40
res.list.Succ.k5.p1.Succ0.40 <- fun.run.sim.Succ(A=20,B=20,Alandscape.LIST.init=Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=5,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ=0.40,dist.vec,array.i,array.j)

pdf("./figs/Succ0.40.Coex.pdf")
par(mfrow=c(2,2),mar=rep(2,4),oma=rep(1,4))
hist((res.list.Succ.k5.p1.Succ0.40[[1]][[20+1]]),main="Early",freq=F)
hist((res.list.Succ.k5.p1.Succ0.40[[2]][[20+1]]),main="Late",freq=F)
image(res.list.Succ.k5.p1.Succ0.40[[1]][[20+1]])
image(res.list.Succ.k5.p1.Succ0.40[[2]][[20+1]])
dev.off()

saveRDS(res.list.Succ.k5.p1.Succ0.40,file="./output/res.list.Succ.k5.p1.Succ0.40.rds")


### plot to explore coexistence
file.name.succ <- c("./output/res.list.Succ.k5.p1.Succ0.95.rds","./output/res.list.Succ.k5.p1.Succ0.85.rds",
  "./output/res.list.Succ.k5.p1.Succ0.5.rds","./output/res.list.Succ.k5.p1.Succ0.40.rds",
  "./output/res.list.Succ.k5.p1.Succ0.35.rds","./output/res.list.Succ.k5.p1.Succ0.30.rds",
  "./output/res.list.Succ.k5.p1.Succ0.25.rds","./output/res.list.Succ.k5.p1.Succ0.20.rds",
  "./output/res.list.Succ.k5.p1.Succ0.15.rds","./output/res.list.Succ.k5.p1.Succ0.1.rds",
                    "./output/res.list.Succ.k5.p1.Succ0.05.rds")

res.coex.succ.list <- lapply(file.name.succ,FUN=readRDS)
vec.param.succ <- c(0.95,0.85,0.5,0.4,0.35,0.3,0.25,0.2,0.15,0.05)

pdf("./figs/Succ.Coexist.pdf")
plot(0,xlim=c(0,1),ylim=c(0,1),pch="",xlab="Succession rate",ylab="Early competitive ability")
for( i in 1:length(res.coex.succ.list)){
res.hist <-  hist(res.coex.succ.list[[i]][[1]][[21]],breaks=(0:10)/10,plot=F)
points(rep(vec.param.succ[i],10),(0:9)/10+0.05,cex=res.hist$count/500)}
abline(v=0.25,col="red")
dev.off()       





###############################################################
###############################################################
#### READ OUTPUT FROM CLUSTER AND PLOT


### NEED to read all outputs from cluster.

## read file in output cluster
list.files(path="./output_cluster",pattern="Succ")

### create vector names of files to read from cluster for the plots
names.p0.disp20.dist0.1 <- c( "res.list.k100.p0.disp20.rds" ,  "res.list.k5.p0.disp20.rds", "res.list.k1.p0.disp20.rds"
                             ,"res.list.k.001.p0.disp20.rds" )
names.p0.5.disp20.dist0.1 <- c( "res.list.k100.p0.5.disp20.rds" ,  "res.list.k5.p0.5.disp20.rds", "res.list.k1.p0.5.disp20.rds"
                               ,"res.list.k.001.p0.5.disp20.rds" )
names.p1.disp20.dist0.1 <- c( "res.list.k100.p1.disp20.rds" ,  "res.list.k5.p1.disp20.rds", "res.list.k1.p1.disp20.rds"
                             ,"res.list.k.001.p1.disp20.rds" )

names.p0.disp20.dist0.5 <- c( "res.list.k100.p0.disp20.dist0.5.rds" ,  "res.list.k5.p0.disp20.dist0.5.rds"
                             , "res.list.k1.p0.disp20.dist0.5.rds" ,"res.list.k.001.p0.disp20.dist0.5.rds" )
names.p0.5.disp20.dist0.5 <- c( "res.list.k100.p0.5.disp20.dist0.5.rds" ,  "res.list.k5.p0.5.disp20.dist0.5.rds"
                               , "res.list.k1.p0.5.disp20.dist0.5.rds" ,"res.list.k.001.p0.5.disp20.dist0.5.rds" )
names.p1.disp20.dist0.5 <- c( "res.list.k100.p1.disp20.dist0.5.rds" ,  "res.list.k5.p1.disp20.dist0.5.rds"
                             , "res.list.k1.p1.disp20.dist0.5.rds" ,"res.list.k.001.p1.disp20.dist0.5.rds" )


names.p0.disp2.dist0.1 <- c( "res.list.k100.p0.disp2.rds" ,  "res.list.k5.p0.disp2.rds", "res.list.k1.p0.disp2.rds"
                            ,"res.list.k.001.p0.disp2.rds" )
names.p0.5.disp2.dist0.1 <- c( "res.list.k100.p0.5.disp2.rds" ,  "res.list.k5.p0.5.disp2.rds", "res.list.k1.p0.5.disp2.rds"
                              ,"res.list.k.001.p0.5.disp2.rds" )
names.p1.disp2.dist0.1 <- c( "res.list.k100.p1.disp2.rds" ,  "res.list.k5.p1.disp2.rds", "res.list.k1.p1.disp2.rds"
                            ,"res.list.k.001.p1.disp2.rds" )

names.p0.disp.2.dist0.1 <- c( "res.list.k100.p0.disp.2.rds" ,  "res.list.k5.p0.disp.2.rds", "res.list.k1.p0.disp.2.rds"
                             ,"res.list.k.001.p0.disp.2.rds" )
names.p0.5.disp.2.dist0.1 <- c( "res.list.k100.p0.5.disp.2.rds" ,  "res.list.k5.p0.5.disp.2.rds", "res.list.k1.p0.5.disp.2.rds"
                               ,"res.list.k.001.p0.5.disp.2.rds" )
names.p1.disp.2.dist0.1 <- c( "res.list.k100.p1.disp.2.rds" ,  "res.list.k5.p1.disp.2.rds", "res.list.k1.p1.disp.2.rds"
                             ,"res.list.k.001.p1.disp.2.rds" )



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

###########################
##### SUCCESSION
names.p0.disp2.dist0.1.Succ0.25 <- c(  "res.list.Succ.k100.p0.Succ0.25.rds" ,
                                     "res.list.Succ.k5.p0.Succ0.25.rds"
                                     , "res.list.Succ.k1.p0.Succ0.25.rds" ,
                                     "res.list.Succ.k.001.p0.Succ0.25.rds")
names.p0.disp2.dist0.1.Succ0.5 <- c(  "res.list.Succ.k100.p0.Succ0.5.rds" ,
                                   "res.list.Succ.k5.p0.Succ0.5.rds"
                                     , "res.list.Succ.k1.p0.Succ0.5.rds" ,
                                   "res.list.Succ.k.001.p0.Succ0.5.rds")
names.p0.disp2.dist0.1.Succ0.85 <- c(  "res.list.Succ.k100.p0.Succ0.85.rds" ,
                                   "res.list.Succ.k5.p0.Succ0.85.rds"
                                     , "res.list.Succ.k1.p0.Succ0.85.rds" ,
                                   "res.list.Succ.k.001.p0.Succ0.85.rds")

pdf("./figs/Firstsimul.F2.SUCC.pdf")
par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(3,4,2,1))
lapply(names.p0.disp2.dist0.1.Succ0.25,fun.plot.grad.quant.cluster.Succ,
       path="output_cluster",lwd1=0.0)
lapply(names.p0.disp2.dist0.1.Succ0.5,fun.plot.grad.quant.cluster.Succ,
       path="output_cluster",lwd1=0.0)
lapply(names.p0.disp2.dist0.1.Succ0.85,fun.plot.grad.quant.cluster.Succ,
       path="output_cluster",lwd1=0.0)
mtext("Succ=0.25",3,adj=0.1, cex=0.9, outer=TRUE)
mtext("Succ=0.5",3,adj=0.5, cex=0.9, outer=TRUE)
mtext("Succ=0.85",3,adj=0.9, cex=0.9, outer=TRUE)
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

pdf("./figs/Firstsimul.image.F.2.pdf")
par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(3,4,2,1))
lapply(names.p0.disp.2.dist0.1 ,fun.image.landscape,
       path="output_cluster")
lapply(names.p0.5.disp.2.dist0.1 ,fun.image.landscape,
       path="output_cluster")
lapply(names.p1.disp.2.dist0.1 ,fun.image.landscape,
       path="output_cluster")
mtext("No Premption (P=0)",3,adj=0.1, cex=0.9, outer=TRUE)
mtext("Half  Premption (P=0.5)",3,adj=0.5, cex=0.9, outer=TRUE)
mtext("Premption (P=1)",3,adj=0.9, cex=0.9, outer=TRUE)
mtext("K=100", 2   , adj=0.92,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=5",  2   , adj=0.65,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=1", 2  , adj=0.37,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=0.001", 2   , adj=0.1,padj=-1.5, cex=1.2, outer=TRUE)
dev.off()

#####
# SUCCESSION IN RGB

pdf("./figs/Firstsimul.image.Succ.F2.pdf")
par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(3,4,2,1),cex.axis=0.1)

lapply(names.p0.disp2.dist0.1.Succ0.25,fun.image.landscape.Succ ,
       path="output_cluster")
lapply(names.p0.disp2.dist0.1.Succ0.5,fun.image.landscape.Succ ,
       path="output_cluster")
lapply(names.p0.disp2.dist0.1.Succ0.85,fun.image.landscape.Succ ,
       path="output_cluster")
mtext("Succ=0.25",3,adj=0.1, cex=0.9, outer=TRUE)
mtext("Succ=0.5",3,adj=0.5, cex=0.9, outer=TRUE)
mtext("Succ=0.85",3,adj=0.9, cex=0.9, outer=TRUE)
mtext("K=100", 2   , adj=0.92,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=5",  2   , adj=0.65,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=1", 2  , adj=0.37,padj=-1.5, cex=1.2, outer=TRUE)
mtext("K=0.001", 2   , adj=0.1,padj=-1.5, cex=1.2, outer=TRUE)
mtext("Climatic gradient", 1   , adj=0.5,padj=1.3, cex=1.4, outer=TRUE)
dev.off()

res.name <- names.p0.disp2.dist0.1.Succ0.85[1]
res.list[[1]][[21]]

#### plot abundance
gray.col.vec <- rev(gray.colors(n=length(res.list)))
abun.temp <- fun.gradient.abundance.levels(res.list,t=1,imax=300)
plot(abun.temp[1,],type="l",ylim=c(0,1000),col=gray.col.vec[1])

for (i in 2:length(res.list)){
abun.temp <- fun.gradient.abundance.levels(res.list,t=i,imax=300)
lines(abun.temp[1,],col=gray.col.vec[i])
}


