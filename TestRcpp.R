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
                     0.1, 0.3, K = 10, n = 500)
l_res_K5<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.3, K = 5, n = 500)

l_res_K100_s0.1<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.1, K = 100, n = 500)
l_res_K10_s0.1<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.1, K = 10, n = 500)
l_res_K5_s0.1<- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                     list_init$c_e, list_init$c_l, list_init$c_s,
                     list_init$ss,
                     0.1, 0.1, K = 5, n = 500)


par(mfrow = c(2, 1))
image_landscape(list_init$mat_sp, list_init$c_e, list_init$c_l, list_init$c_s)
image_landscape(l_res$sp, list_init$c_e, list_init$c_l, list_init$c_s)

### NEED TO EVALUATE TIME TO metaeauilibrium

