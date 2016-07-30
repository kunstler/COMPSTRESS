#####################################
#####################################
#### FUNCTIONS TO PREDICT ABUNDANCE FOR A TWO SPECIES
#### CS trade off along teh gradient based on Case et al. 2005

fun.abund.sup <- function(grad.clim,cs,d,f) {
    apply(cbind(rep(0,length(grad.clim)),
                1-apply(cbind(cs*grad.clim+d,rep(1,length(grad.clim))),
                        MARGIN=1,FUN=min)/f),
          MARGIN=1,FUN=max)
}
fun.abund.inf <- function(grad.clim,ci,cs,d,f) {
    apply(cbind(rep(0,length(grad.clim)),
                2/f*apply(cbind(cs*grad.clim+d,rep(1,length(grad.clim))),
                          MARGIN=1,FUN=min) -
                apply(cbind(ci*grad.clim+d,rep(1,length(grad.clim))),
                      MARGIN=1,FUN=min)/f-1),MARGIN=1,FUN=max)
}




####################################################################
####################################################################
####################################################################
###### CELLULAR SIMULATION MODEL FUNCTIONS

require(raster)
### function to deal with reflexion and semi-torus landscape
## semi torus : if i = 1 then i-k = imaxk+1 / if i = imax then i+k = k
## reflexionprint : if i= 1 then i-k =i+k and if i=Nlandscape *imaxi+k = i-k

# torus function return the N neasrest cells with a torus for one dimension i
function.torus <-  function(i, imax, N = 1){
   if((i+N)>imax){it <- c(i-(N:1),i+(1:N))
      it[it>imax] <-  1:(i+N-imax)
      return(it)}
   else{
       if((i-N)<1){
      it <- c(i-(N:1),i+(1:N))
      it[it<1] <-  (imax+(i-N)):imax
      return(it)}
        else {
         it <- c(i-(N:1),i+(1:N))
         return(it)}
       }
   }


## function.torus(i=5,N=1,imax=256)

## return the N nearest cells with a reflexion
function.reflex <-  function(i, imax, N = 1){
   if((i+N)>imax){it <- c(i-(N:1),i+(1:N))
       it[it>imax] <-  (imax):(2*imax-(i+N-1))
       return(it)} else{if((i-N)<1){
                   it <- c(i-(N:1),i+(1:N))
                   it[it<1] <-  (1:(N-i+1))
                   return(it)} else {it <- c(i-(N:1),i+(1:N))
                               return(it)}
                                    }
   }

## function.reflex(i=250,N=1,imax=256)

## FUNCTION CREATING AN ARRAY OF NEIGHBORHOOD CELLS
function.mat.neigcells <-  function(NN,Nlandscape, N = 1){
mat.index.neigh.cells <- matrix(NA,NN*NN*Nlandscape,8)
mat.cell.num <- matrix(1:(NN*NN*Nlandscape), NN, NN*Nlandscape)

for (i in 1:NN)
{
 for (j in 1:(NN*Nlandscape))
   {
n.i <-  function.torus(i,imax=NN)
n.j <- function.reflex(i=j,imax=NN*Nlandscape)
vec.i<-  c(rep(c(n.i[1:N],i,n.i[(N+1):2*N]),N),
           n.i,
           rep(c(n.i[1:N],i,n.i[(N+1):2*N]),N))
vec.j<- c(rep(n.j[1:N],each=2*N+1),
          rep(j,each=2*N),
          rep(n.j[(N+1):(2*N)],each=2*N+1))
mat.index.neigh.cells[mat.cell.num[i,j], ] <- mat.cell.num[cbind(vec.i, vec.j)]
   }
}
return(mat.index.neigh.cells)
}


function.array.neigcells <-  function(NN,Nlandscape){
array.i <- matrix(NA,dim=c(NN,NN*Nlandscape,8))
array.j <- matrix(NA,dim=c(NN,NN*Nlandscape,8))

for (i in 1:NN)
{
 for (j in 1:(NN*Nlandscape))
   {

n.i <-  function.torus(i,N,imax=NN)
n.j <- function.reflex(i=j,N,imax=NN*Nlandscape)
array.i[i,j,]<-  c(rep(c(n.i[1:N],i,n.i[(N+1):2*N]),N),n.i,
                   rep(c(n.i[1:N],i,n.i[(N+1):2*N]),N))
array.j[i,j,]<- c(rep(n.j[1:N],each=2*N+1),rep(j,each=2*N),
                  rep(n.j[(N+1):(2*N)],each=2*N+1))

   }
}
return(list(array.i,array.j))
}


### function to disperse all seed from neigborhood into a given cell
## disp.fun function with both the dispersal model and the fecundity
## then draw in a poisson distribution form average seed

## function to get the matrix value
function.return.cell <- function(n,i,j,mat) return(mat[i[n],j[n]])


### PROBABLY NOT THE MOST EFFEICIENT TO IMPROVE /THERE WAS NEED TO INCLUDE THE INDIVIDUAL PRESENT IN THE CELL Succ
function.disperse.Succ.to.ij <- function(i,j,Alandscape.LIST,array.i,array.j){
ccc.E.vec <-Alandscape.LIST[[1]][cbind(array.i[i,j,], array.j[i,j,])]
ccc.L.vec <-Alandscape.LIST[[2]][cbind(array.i[i,j,], array.j[i,j,])]
return(list(E = as.vector(na.exclude(c(Alandscape.LIST[[1]][i,j],ccc.E.vec))),
            L = as.vector(na.exclude(c(Alandscape.LIST[[2]][i,j],ccc.L.vec)))))
}



#########
#########
## Function to compute proba of establishment depending on the traits

## based on continous law function
function.law <- function(x,y,K){
return(1/(1+exp(-K*(x-y))))
}


### compute proba of winning all pairs  competitive events
function.compet <-  function(vec.t,K){
mt <- outer(vec.t,vec.t,FUN=function.law,K=K)
diag(mt) <-  1
return(apply(mt,MARGIN=1,FUN=prod)/sum(apply(mt,MARGIN=1,FUN=prod)))
}

colonize.E <- function(vec.c.seed, param.K, Alandscape.LIST){
vec.res <- c(NA,NA)
      if(length(vec.c.seed[[1]])>0) {
          winner <- sample(1:length(vec.c.seed[[1]]),
                           size = 1,
                           prob = function.compet(vec.t=vec.c.seed[[1]],
                                                  K=param.K))
          vec.res[1] <- vec.c.seed[[1]][winner]
          vec.res[2] <- vec.c.seed[[2]][winner]
      }
      else {
         vec.res[1] <-(Alandscape.LIST[[1]][i,j])
         vec.res[2] <-(Alandscape.LIST[[2]][i,j])
      }
return(vec.res)
}

colonize.L <- function(vec.c.seed, param.K, Alandscape.LIST){
vec.res <- c(NA,NA)
     if(length(vec.c.seed[[2]])>0) {
         winner <- sample(1:length(vec.c.seed[[2]]),
                          size = 1,
                          prob = function.compet(vec.t = vec.c.seed[[2]],
                                                 K = param.K))
                      vec.res[1] <- vec.c.seed[[1]][winner]
                      vec.res[2] <- vec.c.seed[[2]][winner]
     } else {
         vec.res[1] <-(Alandscape.LIST[[1]][i,j])
         vec.res[2] <-(Alandscape.LIST[[2]][i,j])
     }
return(vec.res)
}

colonize.NO <- function(vec.c.seed, param.K, Alandscape.LIST){
vec.res <- c(NA,NA)
     if(length(vec.c.seed[[1]])>0) {
         winner <- sample(1:length(vec.c.seed[[1]]),
                          size = 1,
                          prob = function.compet(vec.t=vec.c.seed[[1]],
                                                 K=param.K))
         vec.res[1] <- vec.c.seed[[1]][winner]
         vec.res[2] <- vec.c.seed[[2]][winner]
     } else {
         vec.res[1] <-NA
         vec.res[2] <-NA
      }
return(vec.res)
}


############################################################
## function to update a given cell with new colonizer
### THIS NEZ FUNCTION IS CHANGED TO HAVE TO STATE OF THE CELL EITHER EARLY SUCC OR LATE SUCC
function.colonize.cell.Succ <-  function(i,j,param.K,
                                         Alandscape.LIST,array.i,array.j){
vec.c.seed <- function.disperse.Succ.to.ij(i, j,
                                           Alandscape.LIST,
                                           array.i, array.j)
vec.res <- switch(Alandscape.Succ[i,j],
                  E = colonize.E(vec.c.seed, param.K, Alandscape.LIST),
                  L = colonize.L(vec.c.seed, param.K, Alandscape.LIST),
                  NO = colonize.NO(vec.c.seed, param.K, Alandscape.LIST))
return(vec.res)
}



#######################################################################
### function to kill a cell with disturbance and climate stress
function.kill <-  function(i,j,Alandscape,param.climate.stress,param.dist,fun.clim.morta){
## global dist

 if(is.na(Alandscape[i,j])){
    return(NA)
 }else{
  if(runif(1)<param.dist) {
      return(NA)
  }else{
  ## climate death
  if(runif(1)>fun.clim.morta(Alandscape[i,j],j,param.climate.stress)){
      return(NA)
   }else{
      return(Alandscape[i,j])
   }
  }
 }
}



## climate mortality function
fun.clim.morta1 <- function(c,j,param.climate.stress){
1-c*climate.grad[j]
}

## function.kill(i=1,j=20,Alandscape=AlandscapeE+AlandscapeL,param.climate.stress=NA,param.dist=0.0,fun.clim.morta=fun.clim.morta1)


####
## FUNCTION FOR SUCCESSION FROM E TO L
fun.Succ <-  function(i,j,Alandscape.LIST,param.Succ){
Succ.temp <- Alandscape.LIST[[3]][i,j]
if(Alandscape.LIST[[3]][i,j]=="E"){
    if(runif(1)<param.Succ){Succ.temp <- "L"}
  }
if(Alandscape.LIST[[3]][i,j]=="NO" & !is.na(Alandscape.LIST[[1]][i,j])){
    if(runif(1)<param.Succ){Succ.temp <- "E"}
  }

return(Succ.temp)
}

############################
############################
###### FUNCTION TO UPDATE THE LANDSCAPE ONE STEP


####
fun.update.landscape.Succ<- function(Alandscape.LIST,fun.clim.morta,
                                     param.K,param.climate.stress,
                                     param.dist,param.Succ,array.i,array.j){
for (i in 1:nrow(Alandscape.LIST[[1]])){
   for (j in 1:ncol(Alandscape.LIST[[1]])){
       vec.EL <- function.colonize.cell.Succ(i,j,param.K,
                                             Alandscape.LIST=Alandscape.LIST,
                                             array.i,array.j)
       Alandscape.LIST[[1]][i,j] <- vec.EL[1]
       Alandscape.LIST[[2]][i,j] <- vec.EL[2]

       Alandscape.LIST[[3]][i,j] <- fun.Succ(i,j,Alandscape.LIST,param.Succ)

       morta <- function.kill(i,j,
                              Alandscape= Alandscape.LIST[[1]]+
                                          Alandscape.LIST[[2]],
                              param.climate.stress,
                              param.dist,
                              fun.clim.morta)
       if(is.na(morta)){
        Alandscape.LIST[[1]][i,j] <- NA
        Alandscape.LIST[[2]][i,j] <- NA
        Alandscape.LIST[[3]][i,j] <- "NO"
        }
     }
  }
return(Alandscape.LIST)
}


#########
## function run simulation
fun.run.sim.Succ <- function(A,B,Alandscape.LIST.init,
                             fun.clim.morta=fun.clim.morta1,
                             param.K=1,param.climate.stress=NA,param.dist=0.1,
                             param.Succ,array.i,array.j){

res.listE.temp <- list()
res.listL.temp <- list()
res.list.Succ.temp <- list()
res.listE.temp[[1]] <- Alandscape.LIST.init[[1]]
res.listL.temp[[1]] <- Alandscape.LIST.init[[2]]
res.list.Succ.temp[[1]] <- Alandscape.LIST.init[[3]]
Alandscape.LIST <-  Alandscape.LIST.init

for(j in 1:A){
  for (i in 1:B){
Alandscape.LIST <-  fun.update.landscape.Succ(Alandscape.LIST,fun.clim.morta=fun.clim.morta,
                                              param.K, param.climate.stress,
                                              param.dist,param.Succ,array.i,array.j)
  }
res.listE.temp[[j+1]] <- Alandscape.LIST[[1]]
res.listL.temp[[j+1]] <- Alandscape.LIST[[2]]
res.list.Succ.temp[[j+1]] <- Alandscape.LIST[[3]]
gc()
}
return(list(res.listE.temp,res.listL.temp,res.list.Succ.temp))
}


###################################################
###################################################
### FUNCTION TO ANALYSE OUTPUTS

function.table.levels <-  function(v){
   vec_sp <- factor(v, levels = (0:100)/100)
    table(vec_sp)/sum(table(vec_sp))
}


fun.gradient.table.levels <- function(res.list,t,imax){
sp.abun.grad <- matrix(NA,nrow=101,ncol=imax/10)
    for (i in 1:(imax/10)){
       sp.abun.grad[,i] <- function.table.levels(res.list[[t]][,
                            ((i-1)*10+1):(i*10)])
    }
return(sp.abun.grad)

}

fun.gradient.quantile.levels <- function(res.list,t,imax,probs=
  c(0.05,0.95),nwin=20){
sp.abun.grad <- matrix(NA,nrow=length(probs),ncol=imax/nwin)
    for (i in 1:(imax/nwin)){
       sp.abun.grad[,i] <- quantile(res.list[[t]][,((i-1)*nwin+1):(i*nwin)]
                                    ,probs=probs,na.rm=T)
    }
return(sp.abun.grad)
}

fun.gradient.quantile.levels.mean <- function(res.list,t,imax,probs=
  c(0.05,0.95),nwin=20){
sp.abun.grad <- matrix(NA,nrow=length(probs),ncol=imax/nwin)
    for (i in 1:(imax/nwin)){
      temp <- quantile(res.list[[t[1]]][,((i-1)*nwin+1):(i*nwin)]
                                    ,probs=probs,na.rm=T)

      for(j in t){
       temp <-rbind(temp, quantile(res.list[[j]][,((i-1)*nwin+1):(i*nwin)]
                                    ,probs=probs,na.rm=T))
      }
sp.abun.grad[,i] <- apply(temp,MARGIN=2,FUN=mean)
    }
return(sp.abun.grad)
}

fun.gradient.abundance.levels <- function(res.list,t,imax){
sp.abun.grad <- matrix(NA,nrow=1,ncol=imax/10)
    for (i in 1:(imax/10)){
       sp.abun.grad[,i] <- sum(!is.na(res.list[[t]][,((i-1)*10+1):(i*10)]))
    }
return(sp.abun.grad)

}


#### plot quantile
fun.plot.grad.quant.cluster <-  function(res.name,imax=300,path,lwd1=0.5,lwd2=3){
res.list <- readRDS(file=paste("./",path,"/",res.name,sep=""))
gray.col.vec <- rev(gray.colors(n=length(res.list)))
quant.temp <- fun.gradient.quantile.levels(res.list,t=1,imax=imax)
plot(quant.temp[1,],type="l",ylim=c(0,1),col=gray.col.vec[1],
     ylab="Competitive ability traits",xlab="climate gradient",lwd=lwd1)
lines(quant.temp[2,],col=gray.col.vec[1],lwd=lwd1)

for (i in 2:length(res.list)){
 quant.temp <- fun.gradient.quantile.levels(res.list,t=i,imax=300)
lines(quant.temp[1,],col=gray.col.vec[i],lwd=lwd1)
lines(quant.temp[2,],col=gray.col.vec[i],lwd=lwd1)
}
quant.temp <- fun.gradient.quantile.levels.mean(res.list,t=(length(res.list)-5):length(res.list),imax=imax)
lines(quant.temp[1,],col="red",lwd=lwd2)
lines(quant.temp[2,],col="red",lwd=lwd2)


## add abundance for last step
text(27,0.9,round( sum(!is.na(res.list[[length(res.list)]]))/prod(
  dim(res.list[[length(res.list)]])),4))

}


## PLOT SUCC
#### plot quantile
fun.plot.grad.quant.cluster.Succ <-  function(res.name,imax=300,path,lwd1,lwd2=3){
res.list <- readRDS(file=paste("./",path,"/",res.name,sep=""))
red.col.vec <-  (colorRampPalette(c("grey", "red"))( length(res.list[[1]]) ) )
blue.col.vec <-  (colorRampPalette(c("grey", "blue"))( length(res.list[[1]]) ) )

quant.tempE <- fun.gradient.quantile.levels(res.list[[1]],t=1,imax=imax)
quant.tempL <- fun.gradient.quantile.levels(res.list[[2]],t=1,imax=imax)
plot(quant.tempE[1,],type="l",ylim=c(0,1),col=red.col.vec[1],lwd=0,
     ylab="Competitive ability traits",xlab="climate gradient")
## lines(quant.tempE[2,],col=red.col.vec[1],lwd=lwd1)
## lines(quant.tempL[1,],col=blue.col.vec[1],lwd=lwd1)
## lines(quant.tempL[2,],col=blue.col.vec[1],lwd=lwd1)


## for (i in 2:length(res.list[[1]])){
## quant.tempE <- fun.gradient.quantile.levels(res.list[[1]],t=i,imax=imax)
## quant.tempL <- fun.gradient.quantile.levels(res.list[[2]],t=i,imax=imax)
## lines(quant.tempE[1,],col=red.col.vec[i],lwd=lwd1)
## lines(quant.tempE[2,],col=red.col.vec[i],lwd=lwd1)
## lines(quant.tempL[1,],col=blue.col.vec[i],lwd=lwd1)
## lines(quant.tempL[2,],col=blue.col.vec[i],lwd=lwd1)

## }
quant.tempE <- fun.gradient.quantile.levels.mean(res.list[[1]],t=(length(res.list[[1]])-5):length(res.list[[1]]),imax=imax)
quant.tempL <- fun.gradient.quantile.levels.mean(res.list[[2]],t=(length(res.list[[1]])-5):length(res.list[[1]]),imax=imax)
lines(quant.tempE[1,],type="l",col='red',lwd=lwd2)
lines(quant.tempE[2,],col="red",lwd=lwd2)
lines(quant.tempL[1,],col='blue',lwd=lwd2)
lines(quant.tempL[2,],col='blue',lwd=lwd2)

}


### plot image of landscape
fun.image.landscape <-  function(res.name,path,t=21){
res.list <- readRDS(file=paste("./",path,"/",res.name,sep=""))
image(x=1:nrow(t(res.list[[t]])),y=1:ncol(t(res.list[[t]])),z=t(res.list[[t]])
   ,asp=1,xlab=NA,ylab="climate gradient")
}

### plot image of landscape Succ
fun.image.landscape.Succ <-  function(res.list,t=21){
rast.temp <- stack(raster(res.list[[1]][[t]] ,
                          xmn=1, xmx=ncol(res.list[[1]][[t]]),
                          ymn=1, ymx=nrow(res.list[[1]][[t]])),
                   raster(res.list[[2]][[t]] ,
                          xmn=1, xmx=ncol(res.list[[1]][[t]]),
                          ymn=1,ymx=nrow(res.list[[1]][[t]])),
                   raster (1-(res.list[[1]][[t]] + res.list[[2]][[t]]) ,
                           xmn=1,xmx=ncol(res.list[[1]][[t]]),
                           ymn=1,ymx=nrow(res.list[[1]][[t]])))
plotRGB(rast.temp,scale=1,axes=T ,asp=1)
 }

### function read rds file from cluster

fun.reads.rds.seq <- function(filename,path) {
  object.name <- sub("\\.rds$", "", basename(filename))
  assign(object.name, readRDS(paste("./",path,"/",filename,sep=""))
         ,.GlobalEnv)
}
