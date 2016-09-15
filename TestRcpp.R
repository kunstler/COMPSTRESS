####################################################################
####################################################################
####################################################################
###### Test CELLULAR SIMULATION MODEL in Rcpp to speed up
## R -d valgrind -e "Rcpp::sourceCpp('R/Rcpp/CellAuto.cpp'); source(file='R/RTheoModelFun.R'); fun_debug();"

## Rcpp::sourceCpp('R/Rcpp/CellAuto.cpp'); source(file="R/RTheoModelFun.R"); fun_debug();
## source functions
library(Rcpp)
sourceCpp(file.path('R','Rcpp', 'CellAuto.cpp'))

source(file="R/RTheoModelFun.R")



### NEED TO EVALUATE TIME TO metaequilibrium

## FOR NO spatial
df_sp <- GenerateRandSp(300)
list_init_conv <- InitLandscape(df_sp, NN = 256, Nlandscape = 1,
                           min_ss = 0, max_ss = 0, per_occup = 0.5)

saveRDS(list_init_conv, file = file.path("output", "list_init_conv.rds"))



list_init_conv <- readRDS( file = file.path("output", "list_init_conv.rds"))
conv_K1 <- eval_converg(list_init_conv, Niter = 40 , Nrun = 100,  p_d = 0.1, p_s = 0.15, K = 1)
saveRDS(conv_K1,file = file.path("output", "conv_K1.rds"))

list_init_conv <- readRDS( file = file.path("output", "list_init_conv.rds"))
conv_K5 <- eval_converg(list_init_conv, Niter = 40 , Nrun = 100, p_d = 0.1, p_s = 0.15, K = 5)
saveRDS(conv_K5,file = file.path("output", "conv_K5.rds"))

list_init_conv <- readRDS( file = file.path("output", "list_init_conv.rds"))
conv_K10 <- eval_converg(list_init_conv, Niter = 40 , Nrun = 100, p_d = 0.1, p_s = 0.15, K = 10)
saveRDS(conv_K10,file = file.path("output", "conv_K10.rds"))

list_init_conv <- readRDS( file = file.path("output", "list_init_conv.rds"))
conv_K100 <- eval_converg(list_init_conv, Niter = 40 , Nrun = 100, p_d = 0.1, p_s = 0.15, K = 100)
saveRDS(conv_K100,file = file.path("output", "conv_K100.rds"))

conv_K1 <- readRDS(file = file.path("output", "conv_K1.rds"))
conv_K5 <- readRDS(file = file.path("output", "conv_K5.rds"))
conv_K10 <- readRDS(file = file.path("output", "conv_K10.rds"))
conv_K100 <- readRDS(file = file.path("output", "conv_K100.rds"))


mean_vec_K1<- sapply(conv_K1, function(x, c) mean(c[x[[1]][x[[1]] != 0]+1]), list_init_conv$c_e)
mean_vec_K5<- sapply(conv_K5, function(x, c) mean(c[x[[1]][x[[1]] != 0]+1]), list_init_conv$c_e)
mean_vec_K10<- sapply(conv_K10, function(x, c) mean(c[x[[1]][x[[1]] != 0]+1]), list_init_conv$c_e)
mean_vec_K100<- sapply(conv_K100, function(x, c) mean(c[x[[1]][x[[1]] != 0]+1]), list_init_conv$c_e)

var_vec_K1<- sapply(conv_K1, function(x, c) var(c[x[[1]][x[[1]] != 0]+1]), list_init_conv$c_e)
var_vec_K5<- sapply(conv_K5, function(x, c) var(c[x[[1]][x[[1]] != 0]+1]), list_init_conv$c_e)
var_vec_K10<- sapply(conv_K10, function(x, c) var(c[x[[1]][x[[1]] != 0]+1]), list_init_conv$c_e)
var_vec_K100<- sapply(conv_K100, function(x, c) var(c[x[[1]][x[[1]] != 0]+1]), list_init_conv$c_e)

x11()
par(mfrow = c(4, 1), mar=c(.1,.1,.1,.1))
plot(mean_vec_K1, type = "l", ylim = c(0, 1))
plot(mean_vec_K5, type = "l", ylim = c(0, 1))
plot(mean_vec_K10, type = "l", ylim = c(0, 1))
plot(mean_vec_K100, type = "l", ylim = c(0, 1))

x11()
par(mfrow = c(4, 1), mar=c(.1,.1,.1,.1))
plot(var_vec_K1, type = "l")
plot(var_vec_K5, type = "l")
plot(var_vec_K10, type = "l")
plot(var_vec_K100, type = "l")


## look at spatial pattern is it strange ?
par(mfrow = c(7,7))
lapply(conv_K1, image_landscape, list_init_conv$c_e, list_init_conv$c_l, list_init_conv$c_s)

par(mfrow = c(7,7))
lapply(conv_K5, image_landscape, list_init_conv$c_e, list_init_conv$c_l, list_init_conv$c_s)

par(mfrow = c(7,7))
lapply(conv_K10, image_landscape, list_init_conv$c_e, list_init_conv$c_l, list_init_conv$c_s)

par(mfrow = c(7,7))
lapply(conv_K100, image_landscape, list_init_conv$c_e, list_init_conv$c_l, list_init_conv$c_s)


## FOR spatial
list_init_conv_s <- InitLandscape(df_sp, NN = 256, Nlandscape = 2,
                           min_ss = 0, max_ss = 2, per_occup = 0.5)
saveRDS(list_init_conv_s, file = file.path("output", "list_init_conv_s.rds"))

conv_K1_s <- eval_converg(list_init_conv_s, Niter = 40 , Nrun = 100,  p_d = 0.1, p_s = 0.15, K = 1)
saveRDS(conv_K1_s,file = file.path("output", "conv_K1_s.rds"))

conv_K5_s<- eval_converg(list_init_conv_s, Niter = 40 , Nrun = 100, p_d = 0.1, p_s = 0.15, K = 5)
saveRDS(conv_K5_s,file = file.path("output", "conv_K5_s.rds"))

conv_K10_s <- eval_converg(list_init_conv_s, Niter = 40 , Nrun = 100, p_d = 0.1, p_s = 0.15, K = 10)
saveRDS(conv_K10_s,file = file.path("output", "conv_K10_s.rds"))

conv_K100_s <- eval_converg(list_init_conv_s, Niter = 40 , Nrun = 100, p_d = 0.1, p_s = 0.15, K = 100)
saveRDS(conv_K100_s,file = file.path("output", "conv_K100_s.rds"))


list_init_conv_s <- readRDS( file = file.path("output", "list_init_conv_s.rds"))

conv_K1_s<- readRDS(file = file.path("output", "conv_K1_s.rds"))
conv_K5_s<- readRDS(file = file.path("output", "conv_K5_s.rds"))
conv_K10_s<- readRDS(file = file.path("output", "conv_K10_s.rds"))
conv_K100_s<- readRDS(file = file.path("output", "conv_K100_s.rds"))


mean_vec_K1_s<- sapply(conv_K1_s, function(x, c) mean(c[x[[1]][x[[1]] != 0]+1]), list_init_conv_s$c_e)
mean_vec_K5_s<- sapply(conv_K5_s, function(x, c) mean(c[x[[1]][x[[1]] != 0]]+1), list_init_conv_s$c_e)
mean_vec_K10_s<- sapply(conv_K10_s, function(x, c) mean(c[x[[1]][x[[1]] != 0]+1]), list_init_conv_s$c_e)
mean_vec_K100_s<- sapply(conv_K100_s, function(x, c) mean(c[x[[1]][x[[1]] != 0]+1]), list_init_conv_s$c_e)

var_vec_K1_s<- sapply(conv_K1_s, function(x, c) var(c[x[[1]][x[[1]] != 0]+1]), list_init_conv_s$c_e)
var_vec_K5_s<- sapply(conv_K5_s, function(x, c) var(c[x[[1]][x[[1]] != 0]]+1), list_init_conv_s$c_e)
var_vec_K10_s<- sapply(conv_K10_s, function(x, c) var(c[x[[1]][x[[1]] != 0]+1]), list_init_conv_s$c_e)
var_vec_K100_s<- sapply(conv_K100_s, function(x, c) var(c[x[[1]][x[[1]] != 0]+1]), list_init_conv_s$c_e)

x11()
par(mfrow = c(4, 1), mar=c(.1,.1,.1,.1))
plot(mean_vec_K1_s, type = "l", ylim = c(0, 1))
plot(mean_vec_K5_s, type = "l", ylim = c(0, 1))
plot(mean_vec_K10_s, type = "l", ylim = c(0, 1))
plot(mean_vec_K100_s, type = "l", ylim = c(0, 1))

x11()
par(mfrow = c(4, 1), mar=c(.1,.1,.1,.1))
plot(var_vec_K1_s, type = "l")
plot(var_vec_K5_s, type = "l")
plot(var_vec_K10_s, type = "l")
plot(var_vec_K100_s, type = "l")


## look at spatial pattern is it strange ?
par(mfrow = c(7,7))
lapply(conv_K1_s, image_landscape, list_init_conv_s$c_e, list_init_conv_s$c_l, list_init_conv_s$c_s)

par(mfrow = c(7,7))
lapply(conv_K5_s, image_landscape, list_init_conv_s$c_e, list_init_conv_s$c_l, list_init_conv_s$c_s)

par(mfrow = c(7,7))
lapply(conv_K10_s, image_landscape, list_init_conv_s$c_e, list_init_conv_s$c_l, list_init_conv_s$c_s)

par(mfrow = c(7,7))
lapply(conv_K100_s, image_landscape, list_init_conv_s$c_e, list_init_conv_s$c_l, list_init_conv_s$c_s)



## explore coexistence

df_sp <- GenerateRandSp(300)
list_init_coex<- InitLandscape(df_sp, NN = 256, Nlandscape = 1,
                           min_ss = 0, max_ss = 0, per_occup = 0.5)
saveRDS(list_init_coex, file = file.path("output", "list_init_coex.rds"))

list_coex <- list()
length_p_seq <-  20
p_seq<- seq(0.001, 0.5, length =length_p_seq)
vec_K <-  c(1, 5, 10, 100)
for ( K in 1:length(vec_K)){
 list_coex_K<- list()
   for (d in 1:length_p_seq){
       list_coex_d<- list()
       for (s in 1:length_p_seq){
           res <- UpdateIterR(list_init_coex$mat_sp, list_init_coex$mat_suc,
                            list_init_coex$c_e, list_init_coex$c_l,
                            list_init_coex$c_s,
                            list_init_coex$ss,
                            p_seq[d], p_seq[s], vec_K[K], n = 4000)
           list_coex_d[[s]] <- res
           print("done s")
       }
       list_coex_K[[d]] <- list_coex_d
       print("done d")
   }
 list_coex[[K]] <- list_coex_K
 print("done K")
}
saveRDS(list_coex, file = file.path("output", "list_coex.rds"))

# read output
list_init_coex <- readRDS(file = file.path("output", "list_init_coex.rds"))
list_coex <- readRDS(file = file.path("output", "list_coex.rds"))



par(mfrow = c(1, 2))
image_landscape_e(list_coex[[4]][[2]][[5]], list_init_coex$c_e)
vec_c_e <- list_init_coex$c_e[list_coex[[4]][[2]][[5]][[1]]]
hist(vec_c_e[vec_c_e != -100], breaks = seq(0, 1, length.out = 10))

## Process results
list_coex_sp_n<- list()
list_coex_var_c_e<- list()
for ( K in 1:4){
    mat_n_sp<- matrix(NA, length_p_seq, length_p_seq)
    mat_var_c_e<- matrix(NA, length_p_seq, length_p_seq)
   for (d in 1:length_p_seq){
       for (ss in 1:length_p_seq){
           mat_n_sp[d, ss] <- length(table(list_coex[[K]][[d]][[ss]][[1]]))
           vec_c_e <- list_init_coex$c_e[list_coex[[K]][[d]][[ss]][[1]]]
           mat_var_c_e[d, ss] <- var(vec_c_e[vec_c_e != -100])
       }
   }
list_coex_sp_n[[K]] <- mat_n_sp
list_coex_var_c_e[[K]] <- mat_var_c_e
}

# row distiurbance col succ
pdf(file.path('figs', 'coex_CA_K100.pdf'), width = 14, height = 8)
par(mfrow = c(1, 2))
image(list_coex_sp_n[[4]], col = rev(heat.colors(18)),
                 ylab = 'succession rate', xlab = 'mortality rate')
abline(a = 0, b = 1, col = 'green')
image(list_coex_var_c_e[[4]], col = rev(heat.colors(18)),
                 ylab = 'succession rate', xlab = 'mortality rate')
abline(a = 0, b = 1, col = 'green')
dev.off()

pdf(file.path('figs', 'coex_CA_K10.pdf'), width = 14, height = 8)
par(mfrow = c(1, 2))
image(list_coex_sp_n[[3]], col = rev(heat.colors(18)),
                 ylab = 'succession rate', xlab = 'mortality rate')
abline(a = 0, b = 1, col = 'green')
image(list_coex_var_c_e[[3]], col = rev(heat.colors(18)),
                 ylab = 'succession rate', xlab = 'mortality rate')
abline(a = 0, b = 1, col = 'green')
dev.off()

pdf(file.path('figs', 'coex_CA_K5.pdf'), width = 14, height = 8)
par(mfrow = c(1, 2))
image(list_coex_sp_n[[2]], col = rev(heat.colors(18)),
                 ylab = 'succession rate', xlab = 'mortality rate')
abline(a = 0, b = 1, col = 'green')
image(list_coex_var_c_e[[2]], col = rev(heat.colors(18)),
                 ylab = 'succession rate', xlab = 'mortality rate')
abline(a = 0, b = 1, col = 'green')
dev.off()

pdf(file.path('figs', 'coex_CA_K1.pdf'), width = 14, height = 8)
par(mfrow = c(1, 2))
image(list_coex_sp_n[[1]], col = rev(heat.colors(18)),
                 ylab = 'succession rate', xlab = 'mortality rate')
abline(a = 0, b = 1, col = 'green')
image(list_coex_var_c_e[[1]], col = rev(heat.colors(18)),
                 ylab = 'succession rate', xlab = 'mortality rate')
abline(a = 0, b = 1, col = 'green')
dev.off()

### WHY COEXISTENCE IS IT DIFFERENT FROM THE NON SPATIAL MODEL ??


## COEXISTENCE SEEMS MAXIMUM WHEN MORTALITY RATE EQUAL SUCCESSION RATE ( or smaller above 0.5)
## but also weird pattern NEED ANALITYCAL SOLUTION (NOT YET but numerical solving of ODE confirm that).

## So Run simulation with mortality rate = 0.1 and succession rate = 0.15 (will give higher abundance of late successional species)

df_sp <- GenerateRandSp(100)
list_init <- InitLandscape(df_sp, NN = 256, Nlandscape = 10,
                           min_ss = 0, max_ss = 1, per_occup = 0.5)
saveRDS(list_init, file = file.path("output", "list_init.rds"))



l_res_K100<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.15, K = 100, n = 4000)
saveRDS(l_res_K100,file = file.path("output", "l_res_K100.rds"))
l_res_K10<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.15, K = 10, n = 4000)
saveRDS(l_res_K10,file = file.path("output", "l_res_K10.rds"))
l_res_K5<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.15, K = 5, n = 4000)
saveRDS(l_res_K5,file = file.path("output", "l_res_K5.rds"))
l_res_K1<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.15, K = 1, n = 4000)
saveRDS(l_res_K1,file = file.path("output", "l_res_K1.rds"))

# load
list_init <- readRDS(file = file.path("output", "list_init.rds"))
l_res_K1 <- readRDS(file = file.path("output", "l_res_K1.rds"))
l_res_K5 <- readRDS(file = file.path("output", "l_res_K5.rds"))
l_res_K10 <- readRDS(file = file.path("output", "l_res_K10.rds"))
l_res_K100 <- readRDS(file = file.path("output", "l_res_K100.rds"))


### plot each simulations
pdf(file.path('figs', 'maps_gardients.pdf'))
par(mfrow = c(4, 1))
     image_landscape(l_res_K1,
                     list_init$c_e,
                     list_init$c_l,
                     list_init$c_s)
     mtext(side = 3, text = 'K = 1', line = -4)
     image_landscape(l_res_K5,
                     list_init$c_e,
                     list_init$c_l,
                     list_init$c_s)
     mtext(side = 3, text = 'K = 5', line = -4)

     image_landscape(l_res_K10,
                     list_init$c_e,
                     list_init$c_l,
                     list_init$c_s)
     mtext(side = 3, text = 'K = 10', line = -4)

     image_landscape(l_res_K100,
                     list_init$c_e,
                     list_init$c_l,
                     list_init$c_s)
     mtext(side = 3, text = 'K = 100', line = -4)
dev.off()


## TODO COMPUTE mean and range var of  C_E C_L S along the gradient and number of species and convex hull occupied
