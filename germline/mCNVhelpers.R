
EstimateModeSimple <- function(x) {
  tmpx = x[which(x > 0.5)]
  if (length(tmpx) < 5) {
    return(0)
  }
  
  density_of_x <-  density(tmpx, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  distToMode <- abs(tmpx - mu)
  threshold = max(quantile(distToMode, 0.2), (1 - sqrt(11/12)) / 2)
  if (length(which(distToMode <= threshold) > 40)) {
    return(median(tmpx[which(distToMode <= threshold)]))
  } else {
    lehmanHodges(tmpx[which(distToMode <= threshold)])
  }
}


createMatrixOfLikeliksCoverage <- function(x) {
  
}