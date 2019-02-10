
EstimateModeComplex <- function(x) {
  tmpx = x[which(x > 0.35)]
  if (length(tmpx) < 10) {
    tmpx = x
  }
  matrOfValume <- matrix(rep(tmpx, length(tmpx)), nrow=length(tmpx), byrow=T)
  matrOfValume = sweep(matrOfValume, 1, tmpx, FUN="+") / 2
  values <- apply(combn(tmpx, 2), 2, mean)
  
  density_of_x <-  density(values, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  mu
}
