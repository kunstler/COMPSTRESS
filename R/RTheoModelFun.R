
####################################################################
####################################################################
####################################################################
###### CELLULAR SIMULATION MODEL FUNCTIONS

### function to deal with reflexion and semi-torus landscape
## semi torus : if i = 1 then i-k = imaxk+1 / if i = imax then i+k = k
## reflexion : if i= 1 then i-k =i+k and if i=Nlandscape *imaxi+k = i-k

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
array.i <- array(NA,dim=c(NN ,NN*Nlandscape,8))
array.j <- array(NA,dim=c(NN,NN*Nlandscape,8))
      
for (i in 1:NN)
{
 for (j in 1:(NN*Nlandscape))
   {

n.i <-  function.torus(i,N,imax=dim(Alandscape)[1])
n.j <- function.reflex(i=j,N,imax=dim(Alandscape)[2])
array.i[i,j,]<-  c(rep(c(n.i[1:N],i,n.i[(N+1):2*N]),N),n.i,rep(c(n.i[1:N],i,n.i[(N+1):2*N]),N))
array.j[i,j,]<- c(rep(n.j[1:N],each=2*N+1),rep(j,each=2*N),rep(n.j[(N+1):(2*N)],each=2*N+1))

}}
return(list(array.i,array.j))
}


### function to disperse all seed from neigborhood into a given cell
## disp._fun function with both the dispersal model and the fecundity
## then draw in a poisson distribution form average seed 

## function to get the matrix value
function.return.cell <- function(n,i,j,mat) return(mat[i[n],j[n]])

### PROBABLY NOT THE MOST EFFEICIENT TO IMPROVE
function.disperse.to.ij <- function(i,j,disp.fun,param.DISP,N,Alandscape,dist.vec,array.i,array.j){
ccc.vec <-sapply(1:8,FUN=function.return.cell,i=array.i[i,j,],j=array.j[i,j,],mat=Alandscape)
disp.seed <-  rpois(rep(1,length=length(ccc.vec)),lambda=disp.fun(dist.vec,ccc.vec,param.DISP,N))
return(as.vector(na.exclude(rep(ccc.vec,disp.seed))))
}

## function.disperse.to.ij(i=20,j=3,disp.fun=disp.unif.fun,param.DISP=9,N,Alandscape,dist.vec,array.i,array.j)

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

## function.kill(i=1,j=20,Alandscape,param.climate.stress=NA,param.dist=0.0,fun.clim.morta=fun.clim.morta1)

############################
############################
###### FUNCTION TO UPDATE THE LANDSCAPE ONE STEP

fun.update.landscape<- function(fun.clim.morta,disp.fun,param.DISP,param.K,param.P,N,param.climate.stress,param.dist,dist.vec,array.i,array.j){
for (i in 1:nrow(Alandscape)){ 
   for (j in 1:ncol(Alandscape)){
       Alandscape[i,j]<<- function.colonize.cell(i,j,disp.fun,param.DISP,param.K,param.P,N,Alandscape=Alandscape,dist.vec,array.i,array.j)  
       Alandscape[i,j] <<- function.kill(i,j,Alandscape,param.climate.stress,param.dist,fun.clim.morta)
     }
  }   
}

