## SOLVING differential system of equations

library(deSolve)

parameters <-  c(c = 1, g = 0.3, t = 0.1)
state <- c(E_E = 0.1, E_L = 0.1, L_E = 0.1, L_L = 0.1)

Succ<-function(t, state, parameters) {
   with(as.list(c(state, parameters)),{
     # rate of change
     dE_E<- c*(1- E_E - E_L - L_E - L_L) * (E_E + E_L) + c*L_E *(E_E + E_L) - g* E_E -t * E_E
     dE_L<- g * E_E - c*E_L* (L_E + L_L)  -t * E_L
     dL_E<- c*(1- E_E - E_L - L_E - L_L) * (L_E + L_L) - c*L_E *(E_E + E_L) - g* L_E -t * L_E
     dL_L<- g * L_E + c*E_L* (L_E + L_L)  -t * L_L

     # return the rate of change
     list(c(dE_E, dE_L, dL_E, dL_L))
   })   # end with(as.list ...
 }


times <- seq(0, 100, by = 0.01)
out <- ode(y = state, times = times, func = Succ, parms = parameters)
matplot(out)
diagnostics(out)


times <- seq(0, 5000, length.out = 2)

gg <- seq(0.001, 0.9, length.out = 500)

mat_res <- matrix(NA, nrow = length(gg), ncol = 2)
for (i in seq_len(length(gg))){
parameters <-  c(c = 1, g = gg[i], t = 0.15)
out <- ode(y = state, times = times, func = Succ, parms = parameters)
mat_res[i, ] <- c(sum(out[2, 2:3]), sum(out[2, 4:5]))
}

par(mfrow = c(5, 1), oma = c(0, 0, 0, 0), mar = c(2, 1.5, 1.5, 1))
mat_res[mat_res <1e-10] <-  0
plot(gg, mat_res[, 1], type = 'l')
lines(gg, mat_res[, 2], col = 'red')
abline(v = 0.15, lty = 2)
abline(v = gg[min(which(mat_res[, 1] ==0))], col = 'green')
abline(v = gg[min(which(mat_res[, 2] ==0))], col = 'green')


mat_res <- matrix(NA, nrow = length(gg), ncol = 2)
for (i in seq_len(length(gg))){
parameters <-  c(c = 1, g = gg[i], t = 0.3)
out <- ode(y = state, times = times, func = Succ, parms = parameters)
mat_res[i, ] <- c(sum(out[2, 2:3]), sum(out[2, 4:5]))
}

par(mar = c(2, 1.5, 1.5, 1))
mat_res[mat_res <1e-10] <-  0
plot(gg, mat_res[, 1], type = 'l')
lines(gg, mat_res[, 2], col = 'red')
abline(v = 0.3, lty = 2)
abline(v = gg[min(which(mat_res[, 1] ==0))], col = 'green')
abline(v = gg[min(which(mat_res[, 2] ==0))], col = 'green')

mat_res <- matrix(NA, nrow = length(gg), ncol = 2)
for (i in seq_len(length(gg))){
parameters <-  c(c = 1, g = gg[i], t = 0.5)
out <- ode(y = state, times = times, func = Succ, parms = parameters)
mat_res[i, ] <- c(sum(out[2, 2:3]), sum(out[2, 4:5]))
}

par(mar = c(2, 1.5, 1.5, 1))
mat_res[mat_res <1e-10] <-  0
plot(gg, mat_res[, 1], type = 'l')
lines(gg, mat_res[, 2], col = 'red')
abline(v = 0.5, lty = 2)
abline(v = gg[min(which(mat_res[, 1] ==0))], col = 'green')
abline(v = gg[min(which(mat_res[, 2] ==0))], col = 'green')

mat_res <- matrix(NA, nrow = length(gg), ncol = 2)
for (i in seq_len(length(gg))){
parameters <-  c(c = 1, g = gg[i], t = 0.7)
out <- ode(y = state, times = times, func = Succ, parms = parameters)
mat_res[i, ] <- c(sum(out[2, 2:3]), sum(out[2, 4:5]))
}

par(mar = c(2, 1.5, 1.5, 1))
mat_res[mat_res <1e-10] <-  0
plot(gg, mat_res[, 1], type = 'l')
lines(gg, mat_res[, 2], col = 'red')
abline(v = 0.7, lty = 2)
abline(v = gg[min(which(mat_res[, 1] ==0))], col = 'green')
abline(v = gg[min(which(mat_res[, 2] ==0))], col = 'green')

mat_res <- matrix(NA, nrow = length(gg), ncol = 2)
for (i in seq_len(length(gg))){
parameters <-  c(c = 1, g = gg[i], t = 0.85)
out <- ode(y = state, times = times, func = Succ, parms = parameters)
mat_res[i, ] <- c(sum(out[2, 2:3]), sum(out[2, 4:5]))
}

par(mar = c(2, 1.5, 1.5, 1))
mat_res[mat_res <1e-10] <-  0
plot(gg, mat_res[, 1], type = 'l')
lines(gg, mat_res[, 2], col = 'red')
abline(v = 0.85, lty = 2)
abline(v = gg[min(which(mat_res[, 1] ==0))], col = 'green')
abline(v = gg[min(which(mat_res[, 2] ==0))], col = 'green')


## area of coexistence in g t space

times <- seq(0, 5000, length.out = 2)

gg <- seq(0.01, 0.99, length.out = 500)

mat_E <- mat_L <- matrix(NA, nrow = length(gg), ncol = length(gg))
for (i in seq_len(length(gg))){
  for (j in seq_len(length(gg))){
    parameters <-  c(c = 1, g = gg[i], t = gg[j])
    out <- ode(y = state, times = times, func = Succ, parms = parameters)
    mat_E[i, j] <- sum(out[2, 2:3])
    mat_L[i, j] <- sum(out[2, 4:5])
  }
}

mat_E[mat_E <1e-5] <-  0
mat_L[mat_L <1e-5] <-  0

mat_E_o <-  mat_E
mat_L_o <-  mat_L

mat_E[mat_E <1e-2] <-  0
mat_L[mat_L <1e-2] <-  0

saveRDS(mat_E, file = file.path('output', 'mat_E.rds'))
saveRDS(mat_L, file = file.path('output', 'mat_L.rds'))

## DO THE PLOT OF COEXISTENCE AREA
mat_E <- readRDS(file = file.path('output', 'mat_E.rds'))
mat_L <- saveRDS(file = file.path('output', 'mat_L.rds'))


library(RColorBrewer)

add.alpha <- function(COLORS, ALPHA){
 if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
 RGB <- col2rgb(COLORS, alpha=TRUE)
 RGB[4,] <- round(RGB[4,]*ALPHA)
 NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
 return(NEW.COLORS)
}


pdf(file.path('figs', 'ode_coexistence.pdf'), width = 14, height = 9)
par(mfrow = c(1, 2))
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(mgp = c(2, 0.6, 0))
image(mat_E, col =  add.alpha(colorRampPalette(brewer.pal(9,'Blues'))(12), 0.5),
      axes = FALSE)
  axis(1, at = seq(from = 0.1, to = 0.9, by = 0.1))
  axis(2)
image(mat_L, col =  add.alpha(colorRampPalette(brewer.pal(9,'Greens'))(12), 0.5), add = TRUE, axes = FALSE)
abline(a = 0, b = 1, col = brewer.pal(8,"Set1")[1])
mat <- mat_L
mat[ , ] <-  -1

mat[mat_E >0 & mat_L > 0] <-  1
contour(mat, levels = c(0), add = TRUE)
text(0.05, 0.95, "(a)")
mat_abs <- abs(mat_E - mat_L)
mat_abs[mat == -1] <- 1.5

image(mat_abs, col =  rev(colorRampPalette(brewer.pal(9,'Reds'))(12)),
      axes = FALSE)
abline(a = 0, b = 1, col = brewer.pal(8,"Set1")[2])
  axis(1, at = seq(from = 0.1, to = 0.9, by = 0.1))
mtext("Succession rate g", side = 1, outer = TRUE, line = 2.2)
mtext("Mortality rate g", side = 2, outer = TRUE, line = 2.2)
text(0.05, 0.95, "(b)")
dev.off()
