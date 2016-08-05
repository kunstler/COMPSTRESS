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




par(mfrow = c(4, 1), mar=c(.1,.1,.1,.1))
image_landscape(list_init$mat_sp, list_init$c_e, list_init$c_l, list_init$c_s)
image_landscape(l_res_K5$sp, list_init$c_e, list_init$c_l, list_init$c_s)
image_landscape(l_res_K10$sp, list_init$c_e, list_init$c_l, list_init$c_s)
image_landscape(l_res_K100$sp, list_init$c_e, list_init$c_l, list_init$c_s)


par(mfrow = c(4, 1), mar=c(.1,.1,.1,.1))
image_landscape(list_init$mat_sp, list_init$c_e, list_init$c_l, list_init$c_s)
image_landscape(l_res_K5_s0.1$sp, list_init$c_e, list_init$c_l, list_init$c_s)
image_landscape(l_res_K10_s0.1$sp, list_init$c_e, list_init$c_l, list_init$c_s)
image_landscape(l_res_K100_s0.1$sp, list_init$c_e, list_init$c_l, list_init$c_s)


### NEED TO EVALUATE TIME TO metaeauilibrium

df_sp <- GenerateRandSp(300)
list_init <- InitLandscape(df_sp, NN = 100, Nlandscape = 2,
                           min_ss = 0, max_ss = 1, per_occup = 0.5)



conv_K5 <- eval_converg(list_int, Niter = 20 , Nrun = 200, p_d = 0.1, p_s = 0.3, K = 5)

conv_K1 <- eval_converg(list_int, Niter = 20 , Nrun = 200, p_d = 0.1, p_s = 0.3, K = 1)

conv_K10 <- eval_converg(list_int, Niter = 20 , Nrun = 200, p_d = 0.1, p_s = 0.3, K = 10)

conv_K100 <- eval_converg(list_int, Niter = 20 , Nrun = 200, p_d = 0.1, p_s = 0.3, K = 10)

mean_vec_K1<- sapply(conv_K1, function(x, c) mean(c[x[[1]][x[[1]] != 0]]), list_init$c_e)
mean_vec_K5<- sapply(conv_K5, function(x, c) mean(c[x[[1]][x[[1]] != 0]]), list_init$c_e)
mean_vec_K10<- sapply(conv_K10, function(x, c) mean(c[x[[1]][x[[1]] != 0]]), list_init$c_e)
mean_vec_K100<- sapply(conv_K100, function(x, c) mean(c[x[[1]][x[[1]] != 0]]), list_init$c_e)

par(mfrow = c(4, 1), mar=c(.1,.1,.1,.1))
plot(mean_vec_K1, type = "l")
plot(mean_vec_K5, type = "l")
plot(mean_vec_K10, type = "l")
plot(mean_vec_K100, type = "l")

## convergence seems ok above 1000 or 2000, use 2000 to be safe

## par(mfrow = c(4, 1), mar=c(.1,.1,.1,.1))
## image_landscape(conv_K10[[10]][[1]], list_init$c_e, list_init$c_l, list_init$c_s)
## image_landscape(conv_K10[[12]][[1]], list_init$c_e, list_init$c_l, list_init$c_s)
## image_landscape(conv_K10[[16]][[1]], list_init$c_e, list_init$c_l, list_init$c_s)
## image_landscape(conv_K10[[21]][[1]], list_init$c_e, list_init$c_l, list_init$c_s)


## explore coexistence

df_sp <- GenerateRandSp(300)
list_init <- InitLandscape(df_sp, NN = 256, Nlandscape = 1,
                           min_ss = 0, max_ss = 0, per_occup = 0.5)
list_coex <- list()
p_seq<- seq(0.001, 0.6, length = 20)
vec_K <-  c(1, 5, 10, 100)
for ( K in 1:4){
 list_coex_K<- list()
   for (d in 1:20){
       list_coex_d<- list()
       for (s in 1:20){
           res <- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                            list_init$c_e, list_init$c_l, list_init$c_s,
                            list_init$ss,
                            p_seq[d], p_seq[s], vec_K[K], n = 2000)
           list_coex_d[[s]] <- res
       }
       list_coex_K[[d]] <- list_coex_d
   }
 list_coex[[K]] <- list_coex_K
}



### create random species with a tradeoff between competition stress tolerance and early and late succ compet
## init landscape
df_sp <- GenerateRandSp(300)
list_init <- InitLandscape(df_sp, NN = 256, Nlandscape = 10,
                           min_ss = 0, max_ss = 1, per_occup = 0.5)
l_res_K100<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.3, K = 100, n = 500)
l_res_K10<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.3, K = 10, n = 2000)
l_res_K5<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.3, K = 5, n = 2000)

l_res_K100_s0.1<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.1, K = 100, n = 2000)
l_res_K10_s0.1<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.1, K = 10, n = 2000)
l_res_K5_s0.1<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.1, K = 5, n = 2000)
