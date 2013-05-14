
####################################################################
####################################################################
####################################################################
###### CELLULAR SIMULATION MODEL FUNCTIONS

### function to deal with reflexion and semi-torus landscape
## semi torus : if i = 1 then i-k = imaxk+1 / if i = imax then i+k = k
## reflexionprint : if i= 1 then i-k =i+k and if i=Nlandscape *imaxi+k = i-k

# torus function return the N neasrest cells with a torus for one dimension i
function.torus <-  function(i,N,imax){
   if((i+N)>imax){it <- c(i-(N:1),i+(1:N))
      it[it>imax] <-  1:(i+N-imax)
      return(it)} else{if((i-N)<1){
                  it <- c(i-(N:1),i+(1:N))
                  it[it<1] <-  (imax+(i-N)):imax
                    return(it)} else {it <- c(i-(N:1),i+(1:N))
                                 return(it)}
                                    }
   }


## function.torus(i=5,N=1,imax=256)

## return the N nearest cells with a reflexion
function.reflex <-  function(i,N,imax){
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

function.array.neigcells <-  function(NN,Nlandscape){
array.i <- array(NA,dim=c(NN,NN*Nlandscape,8))
array.j <- array(NA,dim=c(NN,NN*Nlandscape,8))
      
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
## disp._fun function with both the dispersal model and the fecundity
## then draw in a poisson distribution form average seed 

## function to get the matrix value
function.return.cell <- function(n,i,j,mat) return(mat[i[n],j[n]])

### PROBABLY NOT THE MOST EFFEICIENT TO IMPROVE /THERE WAS NEED TO INCLUDE THE INDIVIDUAL PRESENT IN THE CELL
function.disperse.to.ij <- function(i,j,disp.fun,param.DISP,N,
Alandscape ,dist.vec,array.i,array.j){
ccc.vec <-sapply(1:8,FUN=function.return.cell,i=array.i[i,j,],
                  j=array.j[i,j,],mat=Alandscape)
disp.seed <-  rpois(rep(1,length=length(ccc.vec)),lambda=
                    disp.fun(dist.vec,ccc.vec,param.DISP,N))
return(as.vector(na.exclude(c(Alandscape[i,j],rep(ccc.vec,disp.seed)))))
}

### PROBABLY NOT THE MOST EFFEICIENT TO IMPROVE /THERE WAS NEED TO INCLUDE THE INDIVIDUAL PRESENT IN THE CELL Succ
function.disperse.Succ.to.ij <- function(i,j,disp.fun,param.DISP,N,
Alandscape.LIST,dist.vec,array.i,array.j){

ccc.E.vec <-sapply(1:8,FUN=function.return.cell,i=array.i[i,j,],
                  j=array.j[i,j,],mat=Alandscape.LIST[[1]])
ccc.L.vec <-sapply(1:8,FUN=function.return.cell,i=array.i[i,j,],
                  j=array.j[i,j,],mat=Alandscape.LIST[[2]])

disp.seed <-  rpois(rep(1,length=length(ccc.E.vec)),lambda=
                    disp.fun(dist.vec,ccc.vec=ccc.E.vec,param.DISP,N))

return(list(as.vector(na.exclude(c(Alandscape.LIST[[1]][i,j],rep(ccc.E.vec,disp.seed)))),
        as.vector(na.exclude(c(Alandscape.LIST[[2]][i,j],rep(ccc.L.vec,disp.seed))))))
}


#####
## DISPERSAL FUNCTION IN N CELL UNIFORM
disp.unif.fun <-  function(dist.vec=NA,ccc.vec,param.DISP,N){
return(rep(param.DISP/(length((-N):N)^2),length(ccc.vec)))
}

## disp.unif.fun(dist.vec=NA,ccc.vec=1:8,param.DISP=8,N=1)
## matrix distance 8 cells

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

############################################################
## function to update a given cell with new colonizer
### THIS NEZ FUNCTION IS CHANGED TO HAVE TO STATE OF THE CELL EITHER EARLY SUCC OR LATE SUCC
function.colonize.cell <-  function(i,j,disp.fun,param.DISP,param.K,param.P,N,Alandscape,dist.vec,array.i,array.j){
if(!is.na(Alandscape[i,j])){
  if(runif(1)>param.P){
           vec.c.seed <- function.disperse.to.ij(i,j,disp.fun,param.DISP,N,Alandscape,dist.vec,array.i,array.j)
     if(length(vec.c.seed)>0) {return(vec.c.seed[sample(1:length(vec.c.seed),size=1,prob=
                  function.compet(vec.t=vec.c.seed,K=param.K))])} else { return(Alandscape[i,j])}
     }else{return(Alandscape[i,j])} 
   }else{vec.c.seed <- function.disperse.to.ij(i,j,disp.fun,param.DISP,N,Alandscape,dist.vec,array.i,array.j)
         if(length(vec.c.seed)>0) {return(vec.c.seed[sample(1:length(vec.c.seed),size=1,prob=function.compet(vec.t=vec.c.seed,K=param.K))])} else { return(NA)}
        }   
  }  

############################################################
## function to update a given cell with new colonizer
### THIS NEZ FUNCTION IS CHANGED TO HAVE TO STATE OF THE CELL EITHER EARLY SUCC OR LATE SUCC
function.colonize.cell.Succ <-  function(i,j,disp.fun,param.DISP,param.K,param.P,N,
Alandscape.LIST,dist.vec,array.i,array.j){
vec.res <- c(NA,NA)
AlandscapeE <-  Alandscape.LIST[[1]]
AlandscapeL <-  Alandscape.LIST[[2]]
Alandscape.Succ <-  Alandscape.LIST[[3]]


    if(Alandscape.Succ[i,j]=="E"){
           vec.c.seed <- function.disperse.Succ.to.ij(i,j,disp.fun,param.DISP,N,Alandscape.LIST,dist.vec,array.i,array.j)
     if(length(vec.c.seed[[1]])>0) { winner <- (sample(1:length(vec.c.seed[[1]]),size=1,prob=
                  function.compet(vec.t=vec.c.seed[[1]],K=param.K)))
                      vec.res[1] <- vec.c.seed[[1]][winner]
                      vec.res[2] <- vec.c.seed[[2]][winner]  } else { vec.res[1] <-(Alandscape.LIST[[1]][i,j])
                                                                    vec.res[2] <-(Alandscape.LIST[[2]][i,j])}
     }

  if(Alandscape.Succ[i,j]=="L"){
           vec.c.seed <- function.disperse.Succ.to.ij(i,j,disp.fun,param.DISP,N,Alandscape.LIST,dist.vec,array.i,array.j)
     if(length(vec.c.seed[[2]])>0) {winner <- (sample(1:length(vec.c.seed[[2]]),size=1,prob=
                  function.compet(vec.t=vec.c.seed[[2]],K=param.K)))
                      vec.res[1] <- vec.c.seed[[1]][winner]
                      vec.res[2] <- vec.c.seed[[2]][winner]  } else { vec.res[1] <-(Alandscape.LIST[[1]][i,j])
                                                                    vec.res[2] <-(Alandscape.LIST[[2]][i,j])}
     }

   if(Alandscape.Succ[i,j]=="NO"){
          vec.c.seed <- function.disperse.Succ.to.ij(i,j,disp.fun,param.DISP,N,Alandscape.LIST,dist.vec,array.i,array.j)
     if(length(vec.c.seed[[1]])>0) {winner <- (sample(1:length(vec.c.seed[[1]]),size=1,prob=
                  function.compet(vec.t=vec.c.seed[[1]],K=param.K)))
                      vec.res[1] <- vec.c.seed[[1]][winner]
                      vec.res[2] <- vec.c.seed[[2]][winner]
           } else { vec.res[1] <-NA
                    vec.res[2] <-NA}
 
}

return(vec.res)      }   
  





## function.colonize.cell(i=3,j=2,disp.fun=disp.unif.fun,param.DISP=9,param.K=10,param.P=0.5,N=1,Alandscape,dist.vec,array.i,array.j)

#######################################################################
### function to kill a cell with disturbance and climate stress
function.kill <-  function(i,j,Alandscape,param.climate.stress,param.dist,fun.clim.morta){
## global dist

if (is.na(Alandscape[i,j])){return(NA)}else{  if(runif(1)<param.dist) {return(NA)}else{
  ## climate death
  if(runif(1)>fun.clim.morta(Alandscape[i,j],j,param.climate.stress)){return(NA)}else{return(Alandscape[i,j])}
  }}
}



## climate mortality function
fun.clim.morta1 <- function(c,j,param.climate.stress){
1-c*climate.grad[j]
}

## function.kill(i=1,j=20,Alandscape=AlandscapeE+AlandscapeL,param.climate.stress=NA,param.dist=0.0,fun.clim.morta=fun.clim.morta1)


####
## FUNCTION FOR SUCCESSION FROM E TO L
fun.Succ <-  function(i,j,Alandscape.Succ,param.Succ){
Succ.temp <- Alandscape.Succ[i,j]
if(Alandscape.Succ[i,j]=="E"){
    if(runif(1)<param.Succ){Succ.temp <- "L"}
  }
return(Succ.temp)
}

############################
############################
###### FUNCTION TO UPDATE THE LANDSCAPE ONE STEP

fun.update.landscape<- function(Alandscape,fun.clim.morta,disp.fun,param.DISP,param.K,param.P,N,param.climate.stress,
param.dist,dist.vec,array.i,array.j){
for (i in 1:nrow(Alandscape)){ 
   for (j in 1:ncol(Alandscape)){
       Alandscape[i,j]<- function.colonize.cell(i,j,disp.fun,param.DISP,param.K,param.P,N,Alandscape=Alandscape,dist.vec,array.i,array.j)  
       Alandscape[i,j] <- function.kill(i,j,Alandscape,param.climate.stress,param.dist,fun.clim.morta)
     }
  }
return(Alandscape)
}

####
fun.update.landscape.Succ<- function(Alandscape.LIST,fun.clim.morta,disp.fun,param.DISP,
param.K,param.P,N,param.climate.stress,param.dist,param.Succ,dist.vec,array.i,array.j){
for (i in 1:nrow(Alandscape.LIST[[1]])){ 
   for (j in 1:ncol(Alandscape.LIST[[1]])){
       vec.EL <- function.colonize.cell.Succ(i,j,disp.fun,param.DISP,param.K,param.P,N,Alandscape.LIST=Alandscape.LIST,dist.vec,array.i,array.j)  
       Alandscape.LIST[[1]][i,j] <- vec.EL[1]
       Alandscape.LIST[[2]][i,j] <- vec.EL[2]
 
       Alandscape.LIST[[3]][i,j] <- fun.Succ(i,j,Alandscape.Succ=Alandscape.LIST[[3]],param.Succ)

       morta <- function.kill(i,j,Alandscape=Alandscape.LIST[[1]]+Alandscape.LIST[[2]],param.climate.stress,param.dist,fun.clim.morta)
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
fun.run.sim <- function(A,B,Alandscape.init,fun.clim.morta=fun.clim.morta1,disp.fun=disp.unif.fun,param.DISP=2,param.K=1,
                        param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,dist.vec,array.i,array.j){

res.list.temp <- list()
res.list.temp[[1]] <- Alandscape.init
Alandscape <-  Alandscape.init

for(j in 1:A){
  for (i in 1:B){
Alandscape <-  fun.update.landscape(Alandscape=Alandscape,
fun.clim.morta=fun.clim.morta,disp.fun=disp.fun,param.DISP=
param.DISP,param.K=param.K,param.P=param.P,N=N,
  param.climate.stress=param.climate.stress,
  param.dist=param.dist,dist.vec=dist.vec,array.i=array.i,
                                    array.j=array.j)
  }  
res.list.temp[[j+1]] <- Alandscape
gc()
}
return(res.list.temp)
}


#########
## function run simulation
fun.run.sim.Succ <- function(A,B,Alandscape.LIST.init,fun.clim.morta=fun.clim.morta1,
disp.fun=disp.unif.fun,param.DISP=2,param.K=1,param.P=1,N=1,param.climate.stress=NA,param.dist=0.1,param.Succ,dist.vec,array.i,array.j){

res.listE.temp <- list()
res.listL.temp <- list()
res.list.Succ.temp <- list()
res.listE.temp[[1]] <- Alandscape.LIST.init[[1]]
res.listL.temp[[1]] <- Alandscape.LIST.init[[2]]
res.list.Succ.temp[[1]] <- Alandscape.LIST.init[[3]]
Alandscape.LIST <-  Alandscape.LIST.init

for(j in 1:A){
  for (i in 1:B){
Alandscape.LIST <-  fun.update.landscape.Succ(Alandscape.LIST,
fun.clim.morta=fun.clim.morta,disp.fun=disp.fun,param.DISP,param.K,param.P,N,
  param.climate.stress,
  param.dist,param.Succ,dist.vec,array.i,array.j)
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

function.table.levels <-  function(v) table(factor(v,levels=
                                          (0:100)/100))

fun.gradient.table.levels <- function(res.list,t,imax){
sp.abun.grad <- matrix(NA,nrow=101,ncol=imax/10)
    for (i in 1:(imax/10)){
       sp.abun.grad[,i] <- function.table.levels(res.list[[t]][,
                            ((i-1)*10+1):(i*10)])
    }
return(sp.abun.grad)

} 

fun.gradient.quantile.levels <- function(res.list,t,imax,probs=
  c(0.05,0.95)){
sp.abun.grad <- matrix(NA,nrow=length(probs),ncol=imax/10)
    for (i in 1:(imax/10)){
       sp.abun.grad[,i] <- quantile(res.list[[t]][,((i-1)*10+1):(i*10)]
                                    ,probs=probs,na.rm=T)
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
fun.plot.grad.quant.cluster <-  function(res.name,imax=300,path){
res.list <- readRDS(file=paste("./",path,"/",res.name,sep=""))
gray.col.vec <- rev(gray.colors(n=length(res.list)))
quant.temp <- fun.gradient.quantile.levels(res.list,t=1,imax=imax)
plot(quant.temp[1,],type="l",ylim=c(0,1),col=gray.col.vec[1],
     ylab="Competitive ability traits",xlab="climate gradient")
lines(quant.temp[2,],col=gray.col.vec[1])

for (i in 2:length(res.list)){
 quant.temp <- fun.gradient.quantile.levels(res.list,t=i,imax=300)
lines(quant.temp[1,],col=gray.col.vec[i])
lines(quant.temp[2,],col=gray.col.vec[i])
}
lines(quant.temp[1,],col="red")
lines(quant.temp[2,],col="red")
## add abundance for last step
text(27,0.9,round( sum(!is.na(res.list[[length(res.list)]]))/prod(
  dim(res.list[[length(res.list)]])),4))

}


## PLOT SUCC
#### plot quantile
fun.plot.grad.quant.cluster.Succ <-  function(res.name,imax=300,path){
res.list <- readRDS(file=paste("./",path,"/",res.name,sep=""))
red.col.vec <-  (colorRampPalette(c("grey", "red"))( length(res.list[[1]]) ) )
blue.col.vec <-  (colorRampPalette(c("grey", "blue"))( length(res.list[[1]]) ) )

quant.tempE <- fun.gradient.quantile.levels(res.list[[1]],t=1,imax=imax)
quant.tempL <- fun.gradient.quantile.levels(res.list[[2]],t=1,imax=imax)
plot(quant.tempE[1,],type="l",ylim=c(0,1),col=red.col.vec[1],
     ylab="Competitive ability traits",xlab="climate gradient")
lines(quant.tempE[2,],col=red.col.vec[1])
lines(quant.tempL[1,],col=blue.col.vec[1])
lines(quant.tempL[2,],col=blue.col.vec[1])


for (i in 2:length(res.list[[1]])){
quant.tempE <- fun.gradient.quantile.levels(res.list[[1]],t=i,imax=imax)
quant.tempL <- fun.gradient.quantile.levels(res.list[[2]],t=i,imax=imax)
lines(quant.tempE[1,],col=red.col.vec[i])
lines(quant.tempE[2,],col=red.col.vec[i])
lines(quant.tempL[1,],col=blue.col.vec[i])
lines(quant.tempL[2,],col=blue.col.vec[i])

}
## add abundance for last step

}


### plot image of landscape
fun.image.landscape <-  function(res.name,path,t=length(res.list)){
res.list <- readRDS(file=paste("./",path,"/",res.name,sep=""))
image(x=1:nrow(res.list[[t]]),y=1:ncol(res.list[[t]]),z=res.list[[t]]
   ,asp=1,xlab=NA,ylab="climate gradient")
}

### function read rds file from cluster

fun.reads.rds.seq <- function(filename,path) {
  object.name <- sub("\\.rds$", "", basename(filename))
  assign(object.name, readRDS(paste("./",path,"/",filename,sep=""))
         ,.GlobalEnv)
}
