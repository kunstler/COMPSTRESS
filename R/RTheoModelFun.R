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



### PROBABLY NOT THE MOST EFFEICIENT TO IMPROVE /THERE WAS NEED TO INCLUDE THE INDIVIDUAL PRESENT IN THE CELL Succ
return.val.neigh <- function(select, val, mat.neigh){
matrix(val[as.vector(mat.neigh[select, ])], sum(select) , 8)
}

sample.mat <-  function(x,y){
  res <-  NA
  select <- !is.na(x)
  res <- tryCatch(sample(x[select], size = 1, prob = y[select]),
             error = function(e) NA)
return(res)
}

#########
#########
## Function to compute proba of establishment depending on the traits

## based on continous law function
function.law <- function(x,y,K){
return(1/(1 + exp(-K*(x - y))))
}


### compute proba of winning all pairs  competitive events
function.compet <-  function(vec.t,K){
 require(matrixStats)
  mt <- outer(vec.t, vec.t, FUN = function.law, K = K)
  diag(mt) <-  1
  pp <- rowProds(mt, na.rm = TRUE)
  pp[is.na(vec.t)] <-  NA
  return(pp/sum(pp, na.rm = TRUE))
}

colonize.E <- function(mat.state, Ncells, mat.neigh, df.sp, param.K){
  select <- mat.state[, "Succ"] ==2
  mat.neigh.val <- cbind(mat.state[select, "E"],
                         return.val.neigh(select, mat.state[, "E"], mat.neigh))
  mat.prob <- t(apply(mat.neigh.val, MARGIN = 1, function.compet, K = param.K))
  new.val <- mapply(sample.mat, split(mat.neigh.val, row(mat.neigh.val)),
                                     split(mat.prob, row(mat.prob)))
  mat.state[select, "E"] <- new.val
  mat.state[select, "L"] <- df.sp[match(new.val, df.sp$E), "L"]
  return(mat.state[ , 1:2])
}

colonize.L <- function(mat.state, Ncells, mat.neigh, df.sp, param.K){
  if(sum(mat.state[, "Succ"] ==3) > 0){
     select <- mat.state[, "Succ"] ==3
     mat.neigh.val <- cbind(mat.state[select, "L"],
                            return.val.neigh(select, mat.state[, "L"], mat.neigh))
     mat.prob <- t(apply(mat.neigh.val, MARGIN = 1, function.compet, K = param.K))
     new.val <- mapply(sample.mat, split(mat.neigh.val, row(mat.neigh.val)),
                                        split(mat.prob, row(mat.prob)))
     mat.state[select, "L"] <- new.val
     mat.state[select, "E"] <- df.sp[match(new.val, df.sp$L), "E"]
   }
  return(mat.state[ , 1:2])
}

colonize.NO <- function(mat.state, Ncells, mat.neigh, df.sp, param.K){

   select <- mat.state[, "Succ"] ==1
   mat.sel <- mat.state[select, ]
   mat.neigh.val <- return.val.neigh(select, mat.state[, "E"], mat.neigh)
   mat.prob <- t(apply(mat.neigh.val, MARGIN = 1, function.compet, K = param.K))
   new.val <- mapply(sample.mat, split(mat.neigh.val, row(mat.neigh.val)),
                                      split(mat.prob, row(mat.prob)))
   mat.sel[ , "E"] <- new.val
   mat.sel[!is.na(new.val), "L"] <- df.sp[match( new.val[!is.na(new.val)],
                                                df.sp$E), "L"]
   mat.sel[mat.sel[ , "Succ"] ==1 & !is.na(mat.sel[, "E"]) , "Succ"] <- 2
   mat.state[select, ] <- mat.sel
   return(mat.state[, 1:3])
}


kill.cells <- function(mat.state, Ncells,
                       param.dist){

  select <- (1:Ncells)[mat.state[, "Succ"] !=1]
  mat.sel <- mat.state[select, ]
  prob.dead.stress <- apply(mat.sel, MARGIN = 1, sum)
  dead.dist <- runif(length(select))<param.dist
  dead.stress <- runif(length(select)) > prob.dead.stress
  dead.dist | dead.stress
  mat.sel[dead.dist | dead.stress, 1:2] <-  NA
  mat.sel[dead.dist | dead.stress, 3] <-  1
  mat.state[select, ] <- mat.sel
  return(mat.state)
}


succ.cells <- function(mat.state, Ncells,
                       param.succ){
if(sum(mat.state[, "Succ"] ==2)>0){
  select <- (1:Ncells)[mat.state[, "Succ"] ==2]
  succ.L <- runif(length(select))<param.succ
  mat.state[select, "Succ"][succ.L] <- 3
}
  return(mat.state[, 'Succ'])
}

update.cells <- function(mat.state, mat.neigh, Ncells, df.sp,  param.K, param.dist, param.succ){
  mat.state[ , 1:2] <- colonize.E(mat.state, Ncells, mat.neigh, df.sp, param.K)
  mat.state[ , 1:2] <- colonize.L(mat.state, Ncells, mat.neigh, df.sp, param.K)
  mat.state[ , 1:3] <- colonize.NO(mat.state, Ncells, mat.neigh, df.sp, param.K)
  mat.state[ , "Succ"] <- succ.cells(mat.state, Ncells, param.succ)
  mat.state <- kill.cells(mat.state, Ncells, param.dist)
  return(mat.state)
}

#########
## function run simulation
fun.run.sim.succ <- function(A,B, mat.state, mat.neigh, Ncells, df.sp, param.K, param.dist, param.succ){
res.list.temp <- list()
res.list.temp[[1]] <- mat.state
for(j in 1:A){
  for (i in 1:B){
mat.state <- update.cells(mat.state, mat.neigh, Ncells, df.sp, param.K, param.dist, param.succ)
  }
res.list.temp[[j+1]] <- mat.state
}
return(res.list.temp)
}


generate.rand.sp <- function(nsp, Nval = 10000){
   if(nsp > 5000) stop("error nsp must be smaller than 5000")
   E <- matrix(rep(0:Nval, Nval + 1), Nval +1, Nval+1)
   L <- matrix(rep(0:Nval, Nval + 1), Nval + 1, Nval +1, byrow= TRUE)
   S <- (E+L)
   E.r <-  E[rev(0:(Nval + 1)), ]
   L.r <-  L[rev(0:(Nval + 1)), ]
   S.r <-  S[rev(0:(Nval + 1)), ]
   E.r[upper.tri(E.r,diag = FALSE)] <-  NA
   L.r[upper.tri(E.r,diag = FALSE)] <-  NA
   S.r[upper.tri(E,diag = FALSE)] <-  NA
   E.b <- E.r[rev(0:(Nval + 1)), ]
   L.b <- L.r[rev(0:(Nval + 1)), ]
   S.b <- S.r[rev(0:(Nval + 1)), ]

   df <- na.omit(data.frame(E = as.vector(E.b)/Nval,
                            L = as.vector(L.b)/Nval,
                            S = as.vector(S.b)/Nval))
  df.sp <- df[sample(1:dim(df)[1], nsp+4), ]
  df.sp <- df.sp[!duplicated(df.sp$E) & !duplicated(df.sp$L), ]
   return(df.sp)
}


plot.sp <- function(df){
   par(mfrow=c(2,2))
   x   <- cbind(df$E , df$L)
   data <- as.data.frame(x)
   names(data) <- c("X","Y")
   colors  <- densCols(x)
   plot(x, col=colors, pch=20,cex=0.25,xlab="C.E",ylab="C.L")
   abline(a = 1, b = -1, col = 'red')
   x   <- cbind(1 - df$E - df$L, df$L)
   data <- as.data.frame(x)
   names(data) <- c("X","Y")
   colors  <- densCols(x)
   plot(x, col=colors, pch=20,cex=0.25,xlab="S",ylab="C.L")
   abline(a = 1, b = -1, col = 'red')
   x   <- cbind(df$E, 1 - df$E - df$L)
   data <- as.data.frame(x)
   names(data) <- c("X","Y")
   colors  <- densCols(x)
   plot(x, col=colors, pch=20,cex=0.25,xlab="C.E",ylab="S")
   abline(a = 1, b = -1, col = 'red')
}

###################################################
###################################################
### FUNCTION TO ANALYSE OUTPUTS

function.table.levels <-  function(v, nn = 1){
   vec_sp <- table(cut(v[,nn], breaks = (0:100)/100))
   vec_sp/sum(vec_sp, na.rm = TRUE)
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

