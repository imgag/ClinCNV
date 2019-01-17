
sdsOfRegions <- apply(coverage, 1, sd)
potentiallyPolymorphicRegions <- which(sdsOfRegions > quantile(sdsOfRegions, 0.9) & !bedFile[,1] %in% c("chrX","chrY"))

coverageForClustering = coverage[-potentiallyPolymorphicRegions,]
n = 5
coverageForClustering = (apply(coverageForClustering, 2, function(x) tapply(x, ceiling(seq_along(x) / n), mean)))

corMatrix <- cor(coverageForClustering)

distMatrix <- 1 - corMatrix

for (i in 2:100) {
  hc = hclust(as.dist(distMatrix), method="ward.D")
  memb = cutree(hc, k = i)
  if (min(table(memb)) < 50) {
    break
  }
}