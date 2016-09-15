#####################################
#####################################
## Functions for CellAuto


##########
## Init landscape

InitLandscape <- function(df_sp, NN = 100, Nlandscape = 1, min_ss = 0, max_ss = 0, per_occup = 0.5){
   matrix_sp <- matrix(0, nrow=NN,ncol=Nlandscape*NN)
   matrix_suc <- matrix(1, nrow=NN,ncol=Nlandscape*NN)
   climate_grad <-  seq(from= min_ss,to= max_ss,length=Nlandscape*NN) # not used yet
   ## Select cells to init landscape with species
   init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*per_occup))
   ## fill with species
   nsp <- nrow(df_sp)
   matrix_sp[init.temp] <- sample(1:nsp, length(init.temp), replace = TRUE)

   ## species param
   sp_e <- c(-100, df_sp$E)
   sp_l <- c(-100, df_sp$L)
   sp_s <- c(-100, df_sp$S)

   return(list(mat_sp = matrix_sp, mat_suc = matrix_suc, ss = climate_grad,
               c_e = sp_e, c_l = sp_l, c_s = sp_s))
}

## generate species
GenerateRandSp <- function(nsp, Nval = 1000){
   if(nsp > 400) stop("error nsp must be smaller than 400")
   E <- matrix(rep(0:Nval, Nval + 1), Nval +1, Nval+1)/Nval
   L <- matrix(rep(0:Nval, Nval + 1), Nval + 1, Nval +1, byrow= TRUE)/Nval
   S <- (E+L)
   E.r <-  E[rev(0:(Nval + 1)), ]
   L.r <-  L[rev(0:(Nval + 1)), ]
   S.r <-  S[rev(0:(Nval + 1)), ]
   E.r[upper.tri(E.r,diag = FALSE)] <-  NA
   L.r[upper.tri(L.r,diag = FALSE)] <-  NA
   S.r[upper.tri(S.r,diag = FALSE)] <-  NA
   E.b <- E.r[rev(0:(Nval + 1)), ]
   L.b <- L.r[rev(0:(Nval + 1)), ]
   S.b <- S.r[rev(0:(Nval + 1)), ]

   df <- na.omit(data.frame(E = as.vector(E.b),
                            L = as.vector(L.b),
                            S = as.vector(S.b)))
  df.sp <- df[sample(1:dim(df)[1], nsp+50), ]
  df.sp <- df.sp[!duplicated(df.sp$E) & !duplicated(df.sp$L), ]
if(nrow(df.sp) > nsp) df.sp <- df.sp[sample(1:nrow(df.sp), nsp), ]
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

return_mat_fill_c_raster<- function(m, cc){
  require(raster)
  d <-  dim(m)
  cc[1] <- NA
  mat <- matrix(cc[m+1], d[1], d[2])
  return(raster(mat,
                xmn=1, xmx=d[2],
                ymn=1, ymx=d[1]))
}

return_mat_raster<- function(m){
  require(raster)
  d <-  dim(m)
  return(raster(m,
                xmn=1, xmx=d[2],
                ymn=1, ymx=d[1]))
}

### plot image of landscape Succ

image_return_mat<-  function(list_res, ...){
    rast_e <-  return_mat_raster(list_res)
    image(rast_e,axes=FALSE ,asp=1, ...)
 }

image_landscape_e<-  function(list_res, c_e){
    rast_e <-  return_mat_fill_c_raster(list_res[[1]], c_e)
    image(rast_e,axes=FALSE ,asp=1)
 }

image_landscape <-  function(list_res, c_e, c_l, c_s, ...){
    rast_e <-  return_mat_fill_c_raster(list_res[[1]], c_e)
    rast_l <-  return_mat_fill_c_raster(list_res[[1]], c_l)
    rast_s <-  return_mat_fill_c_raster(list_res[[1]], c_s)
    rast.temp <- stack(rast_e,
                       rast_l,
                       rast_s)
    plotRGB(rast.temp,scale=1,axes=FALSE ,asp=1, ...)
 }


## test convergence

eval_converg<- function(list_init, Niter, Nrun, p_d, p_s, K){

  list_res <- vector('list', Niter+1)
  list_res[[1]] <- list_init
  print(1)
  list_res[[2]] <- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                       list_init$c_e, list_init$c_l, list_init$c_s,
                       list_init$ss,
                       p_d, p_s, K , n = Nrun)
  for (i in 2:Niter){
     print(i)
  list_res[[i + 1]] <- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                       list_init$c_e, list_init$c_l, list_init$c_s,
                       list_init$ss,
                       p_d, p_s, K , n = i*Nrun)
     ## list_res[[i+1]] <- UpdateIterR(list_res[[i]]$sp, list_res[[i]]$suc,
     ##                   list_init$c_e, list_init$c_l, list_init$c_s,
     ##                   list_init$ss,
     ##                   p_d, p_s, K , n = Nrun)
     ## I do not understand why this is not working ....
  }
  return(list_res)
}


## coexistence

eval_coex<- function(p_d, p_s, K, list_int, Nrun){
res <- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                       list_init$c_e, list_init$c_l, list_init$c_s,
                       list_init$ss,
                       p_d, p_s, K , n = Nrun)
return(res[[1]])
}


## function to process output


table_level <- function(x, nsp){
  v <- rep(0, nsp+1)
  names(v) <- 0:nsp
  r <- table(x)
  v[names(r)] <-  r
  return(v)
}

table_level_row <- function(mm, nsp){

t(apply(mm, MARGIN = 1, table_level, nsp))
}


