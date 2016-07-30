
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
## initialization  Succ TODO wrap up in a function
## Landscapoe matrix
Nlandscape <-  1 #dimension (in Nlandscape*NN) of the long side corressponding to the climatic stress gradient
NN <- 25  ## size of landscape (classicaly 256)
matrix.cells.index <- matrix(1:(NN*NN*Nlandscape), nrow=NN,ncol=Nlandscape*NN)
climate.grad <-  seq(from=0,to=0,length=Nlandscape*NN)

matrix.cells.stress <- matrix(rep(climate.grad , NN),
                              nrow=NN, ncol=Nlandscape*NN,
                              byrow = TRUE)

mat.state.init <- matrix(NA, nrow = (NN*NN*Nlandscape), ncol = 4)
colnames(mat.state.init) <- c("E", "L", "Succ", "Stress")
mat.state.init[, "Succ"] <-  1
mat.state.init[, "Stress"] <-  as.vector(matrix.cells.stress)



## init landscape with species
init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*0.5))

### create random species with a tradeoff between competition stress tolerance and eraly and late succ compet
df.sp <- generate.rand.sp(100)

df.samp <- df.sp[sample(1:nrow(df.sp), round(NN*NN*Nlandscape*0.5), replace = TRUE), ]
plot.sp(df.samp)

## fun to create random early succ and late succ competitive ability
mat.state.init[init.temp, "E"] <- df.samp$E
mat.state.init[init.temp, "L"] <- df.samp$L
mat.state.init[init.temp, "Succ"] <- 2



## create a table of 8 neigborhood cells index with torus and reflection
mat.neigh <- function.mat.neigcells(NN,Nlandscape)

Rprof("/tmp/succ.prof")
res <- fun.run.sim.succ(A = 3,B = 30, mat.state = mat.state.init, mat.neigh, Ncells = nrow(mat.neigh), df.sp = df.sp, param.K = 5, param.dist = 0.1, param.succ = 0.25)
Rprof(NULL)
head(summaryRprof("/tmp/succ.prof")$by.self, 10)

system.time(
    res <- fun.run.sim.succ(A = 4,B = 100, mat.state = mat.state.init, mat.neigh, Ncells = nrow(mat.neigh), df.sp = df.sp, param.K = 5, param.dist = 0.1, param.succ = 0.25)
)
hist(res[[5]][, 1])
hist(res[[5]][, 2])
plot(res[[5]][,1], res[[5]][,2])

sapply(res, function(x) table(x[,3]))
table(res[[2]][, 3])

##################################
##################################
#### RUN A SIMULATION SUCCESSION TES  FOR A B time step
## K 100 P 1 Succ 0.25


fun.coex <- function(param.S){
res <- fun.run.sim.succ(A = 1,B = 400, mat.state = mat.state.init, mat.neigh,
                        Ncells = nrow(mat.neigh), df.sp = df.sp, param.K = 5,
                        param.dist = 0.1, param.succ = 0.25)
return(res[[2]])
}

res_list <- lapply(seq(0.1, 0.5, by = 0.05), fun.coex)
image(sapply(res_list,function.table.levels), col = rev(heat.colors(50)))


