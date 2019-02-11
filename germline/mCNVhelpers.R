
EstimateModeSimple <- function(x) {
  tmpx = x[which(x > 0.3)]
  if (length(tmpx) < 10) {
    return(0)
  }
  
  density_of_x <-  density(tmpx, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  distToMode <- abs(tmpx - mu)
  median(tmpx[which(distToMode < quantile(distToMode, 0.2))])
}
