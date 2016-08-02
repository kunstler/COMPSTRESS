####################################################################
####################################################################
####################################################################
###### Test CELLULAR SIMULATION MODEL in Rcpp to speed up

## source functions
library(Rcpp)
sourceCpp(file.path('R','Rcpp', 'CS_Cell.cpp'))

source(file="R/RTheoModelFun.R")


## init landscape to be moved in R functions

Nlandscape <-  1 #dimension (in Nlandscape*NN) of the long side corressponding to the climatic stress gradient
NN <- 256  ## size of landscape (classicaly 256)
matrix_sp <- matrix(0, nrow=NN,ncol=Nlandscape*NN)
matrix_suc <- matrix(1, nrow=NN,ncol=Nlandscape*NN)
climate.grad <-  seq(from=0,to=0,length=Nlandscape*NN) # not used yet


## Select cells to init landscape with species
init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*0.5))

### create random species with a tradeoff between competition stress tolerance and eraly and late succ compet
df.sp <- generate.rand.sp(100)

plot.sp(df.sp)

nsp <- nrow(df.sp)

matrix_sp[init.temp] <- sample(1:nsp, length(init.temp), replace = TRUE)

head(df.sp)
sp_e <- c(-100, df.sp$E)
sp_l <- c(-100, df.sp$L)

system.time(
list.res <- UpdateIterTest(matrix_sp, matrix_suc, sp_e, sp_l, 500)
)

# test object oriented version

sourceCpp(file.path('R','Rcpp', 'CellAuto.cpp'))


l_res <- UpdateIterR(matrix_sp, matrix_suc, 256, 256, sp_e, 0.1, 10)
