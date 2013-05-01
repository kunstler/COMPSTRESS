
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
Nlandscape <-  2 #dimension (in 2*NN) of the long side corressponding to the climatic stress gradient
NN <-  10 ## size of landscape (classicaly 256)
N <-  1 ## number of neigborhood cell 1 = moor neigborhood
Alandscape <-  matrix(NA,nrow=NN,ncol=Nlandscape*NN)
## climate gradient
climate.grad <-  seq(from=0,to=1,length=Nlandscape*NN)

## init landscape with random draw of competition traits
init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*0.5))
vec.land <- as.vector(Alandscape)
vec.land[init.temp] <- sample((0:100)/100,size=round(NN*NN*Nlandscape*0.5),replace=T)
Alandscape <- matrix(vec.land,nrow=NN,ncol=NN*Nlandscape)

## create a table of 8 neigborhood cells

list.temp <- function.array.neigcells(NN,Nlandscape)

array.i <- list.temp[[1]]
array.j <- list.temp[[2]]
rm(list.temp)
array.i[10,20,]



##################################
##################################
#### RUN A SIMULATION FOR A B time step
## system.time(temp <- fun.udpate.landscape(Alandscape,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=8,param.K=10,param.P=0,N=1,param.climate.stress=NA,param.dist=0.5,dist.vec,array.i,array.j))
A <-  4
B <-  10
res.list <- list()
res.list[[1]] <- Alandscape

system.time(
for(j in 1:A){
for (i in 1:B){
fun.update.landscape(fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=2,param.P=0,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j)

}  
res.list[[j+1]] <- Alandscape

gc()
}
)

length(res.list)

saveRDS(res.list,file="./output/res.list.rds")

####read old RDS
res.list <- readRDS("./output/res.list.rds")


## plot landscape at different time step
image(x=1:nrow(res.list[[1]]),y=1:ncol(res.list[[1]]),z=res.list[[1]],xlim=c(1,ncol(res.list[[1]])),ylim=c(1,ncol(res.list[[1]])))
x11(); image(x=1:nrow(res.list[[1]]),y=1:ncol(res.list[[1]]),z=res.list[[2]],xlim=c(1,ncol(res.list[[1]])),ylim=c(1,ncol(res.list[[1]])))
x11(); image(x=1:nrow(res.list[[1]]),y=1:ncol(res.list[[1]]),z=res.list[[3]],xlim=c(1,ncol(res.list[[1]])),ylim=c(1,ncol(res.list[[1]])))
x11(); image(x=1:nrow(res.list[[1]]),y=1:ncol(res.list[[1]]),z=res.list[[5]],xlim=c(1,ncol(res.list[[1]])),ylim=c(1,ncol(res.list[[1]])))





## proceesing output
function.table.levels <-  function(v) sum((0:100)/100*table(factor(v,levels=(0:100)/100)))/length(v[!is.na(v)])
function.table.levels <-  function(v) table(factor(v,levels=(0:100)/100))

for (i in 1:length(res.list))
{apply(res.list[[i]],FUN=function.table.levels,MARGIN=2)
}

(apply(res.list[[3]],FUN=function.table.levels,MARGIN=2))
for(i in 1:4) image(res.list[[i]])

x11(); image(res.list[[10]])
      ?image
dim(res.list[[1]])
image(res.list[[16]])

length(res.list)


?system.time

?sample
?matrix

mtest <-matrix(1:16,nrow=4,ncol=4)
plot(as.vector(row(mtest)),as.vector(col(mtest)))

fun.tes <- function(mtest){if(mtest[x,y]>5){return(10)}else{return(0)}}
apply(X=1:4,Y=1:4,FUN=fun.tes,mtest)

library(simecol)
?lapply
?mapply
?outer
