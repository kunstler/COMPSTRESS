
####################################################################
####################################################################
####################################################################
###### CELLULAR SIMULATION MODEL

## K 100 5 1 0.001
# p 0
# disp 20

## source functions
source(file="R/RTheoModelFun.R")



##################################
## initialization  Succ
## Landscapoe matrix
Nlandscape <-  1 #dimension (in 2*NN) of the long side corressponding to the climatic stress gradient
NN <- 25  ## size of landscape (classicaly 256)
N <-  1 ## number of neigborhood cell 1 = moor neigborhood
AlandscapeE <-  matrix(NA,nrow=NN,ncol=Nlandscape*NN)
AlandscapeL <- AlandscapeE
Alandscape.Succ <- matrix("NO",nrow=NN,ncol=Nlandscape*NN)

## climate gradient
climate.grad <-  seq(from=0,to=0,length=Nlandscape*NN)

## init landscape with random draw of competition traits
init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*0.5))
vec.land <- as.vector(AlandscapeE)
vec.landE <- as.vector(AlandscapeE)
vec.landL <- as.vector(AlandscapeE)
vec.land.Succ <- as.vector(Alandscape.Succ)

### create random species with a tradeoff between competition stress tolerance and eraly and late succ compet
E <- matrix(rep(0:100, 101), 101, 101)
L <- matrix(rep(0:100, 101), 101, 101, byrow= TRUE)
S <- (E+L)
E.r <-  E[rev(0:101), ]
L.r <-  L[rev(0:101), ]
S.r <-  S[rev(0:101), ]
E.r[upper.tri(E.r,diag = FALSE)] <-  NA
L.r[upper.tri(E.r,diag = FALSE)] <-  NA
S.r[upper.tri(E,diag = FALSE)] <-  NA
E.b <- E.r[rev(0:101), ]
L.b <- L.r[rev(0:101), ]
S.b <- S.r[rev(0:101), ]

df <- na.omit(data.frame(E = as.vector(E.b),
                         L = as.vector(L.b),
                         S = as.vector(S.b)))
num.sp.select <- sample(sample(1:dim(df)[1], 100), round(NN*NN*Nlandscape*0.5), replace = TRUE)

## fun to create random early succ and late succ competitive ability
vec.landE[init.temp] <- df$E[num.sp.select]/100
vec.landL[init.temp] <- df$L[num.sp.select]/100
vec.land[init.temp] <- df$S[num.sp.select]/100

par(mfrow=c(2,2))
x   <- cbind(vec.landE[init.temp] ,vec.landL[init.temp])
data <- as.data.frame(x)
names(data) <- c("X","Y")
colors  <- densCols(x)
plot(x, col=colors, pch=20,cex=0.25,xlab="C.E",ylab="C.L")
abline(a = 1, b = -1, col = 'red')
x   <- cbind(1-vec.land[init.temp] ,vec.landL[init.temp])
data <- as.data.frame(x)
names(data) <- c("X","Y")
colors  <- densCols(x)
plot(x, col=colors, pch=20,cex=0.25,xlab="S",ylab="C.L")
abline(a = 1, b = -1, col = 'red')
x   <- cbind(vec.landE[init.temp] ,1-vec.land[init.temp])
data <- as.data.frame(x)
names(data) <- c("X","Y")
colors  <- densCols(x)
plot(x, col=colors, pch=20,cex=0.25,xlab="C.E",ylab="S")
abline(a = 1, b = -1, col = 'red')


vec.land.Succ[init.temp] <-  rep("E",round(NN*NN*Nlandscape*0.5))
AlandscapeE <- matrix(vec.landE,nrow=NN,ncol=NN*Nlandscape)
AlandscapeL <- matrix(vec.landL,nrow=NN,ncol=NN*Nlandscape)
Alandscape.Succ <- matrix(vec.land.Succ,nrow=NN,ncol=NN*Nlandscape)

Alandscape.LIST.init <-  list(AlandscapeE,AlandscapeL,Alandscape.Succ)

## create a table of 8 neigborhood cells

list.temp <- function.array.neigcells(NN,Nlandscape)

image(matrix(1:625 %in% mat.temp[100,],NN,NN*Nlandscape))
mat.temp <- function.mat.neigcells(NN,Nlandscape)
vec.e <- sample((0:100)/100, size = NN*NN*Nlandscape,replace = TRUE)
neigh.val <- cbind(vec.e,
                   matrix(vec.e[as.vector(mat.temp)], NN*NN*Nlandscape, 8))
mat.prob <- t(apply(neigh.val, MARGIN = 1, function.compet, K = 5))

sample.mat <-  function(x,y){
    sample(x, size = 1, prob = y)
}

sampled.cell <- mapply(sample.mat, split(neigh.val, row(neigh.val)),
                                   split(mat.prob, row(mat.prob)))

array.i <- list.temp[[1]]
array.j <- list.temp[[2]]
rm(list.temp)

matrix(1:prod(dim(Alandscape)),
##################################
##################################
#### RUN A SIMULATION SUCCESSION TES  FOR A B time step
## K 100 P 1 Succ 0.25

Rprof("/tmp/succ.prof")
res.test <- fun.run.sim.Succ(A=1,B=200,
                        Alandscape.LIST.init=Alandscape.LIST.init,
                        fun.clim.morta=fun.clim.morta1,
                        param.K=5,
                        param.climate.stress=NA,
                        param.dist=0.1,
                        param.Succ=0.25,
                        array.i = array.i,
                        array.j = array.j)
Rprof(NULL)
head(summaryRprof("/tmp/succ.prof")$by.self, 10)

fun.coex <- function(param.S){
res <- fun.run.sim.Succ(A=1,B=800,
                        Alandscape.LIST.init=Alandscape.LIST.init,
                        fun.clim.morta=fun.clim.morta1,
                        param.K=5,
                        param.climate.stress=NA,
                        param.dist=0.1,
                        param.Succ=param.S,
                        array.i = array.i,
                        array.j = array.j)
return(res[[1]][[2]])
}

res_list <- lapply(seq(0.1, 0.5, by = 0.05), fun.coex)
image(sapply(res_list,function.table.levels), col = rev(heat.colors(50)))

fun.image.landscape.Succ(res.test, 11)

image(sapply(res.test[[1]],function.table.levels), col = rev(heat.colors(50)))
