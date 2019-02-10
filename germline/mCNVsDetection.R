setwd(opt$folderWithScript)
source(paste0(opt$folderWithScript, "/germline/mCNVhelpers.R"))


signalToNoise = mediansAndSds[[2]][which(!bedFile[,1] %in% c("chrX", "chrY")),1] / (mediansAndSds[[2]][which(!bedFile[,1] %in% c("chrX", "chrY")),2])
signalToNoiseSex = mediansAndSds[[2]][which(bedFile[,1] %in% c("chrX", "chrY")),1] / (mediansAndSds[[2]][which(bedFile[,1] %in% c("chrX", "chrY")),2])
signalToNoiseAll = mediansAndSds[[2]][,1] / mediansAndSds[[2]][,2]
potentiallyPolymorphic <- which(mediansAndSds[[2]][,1] > sqrt(3/2) | 
                                  (signalToNoiseAll < quantile(signalToNoise, 0.05) & !bedFile[,1] %in% c("chrX", "chrY")) |
                                  (signalToNoiseAll < quantile(signalToNoiseSex, 0.05) & bedFile[,1] %in% c("chrX", "chrY"))
                                  )

mixOfSamples <- sample(1:ncol(coverage))
mediansOfPolymorphic = foreach (i=1:nrow(coverage), .combine="c") %dopar% {
  modeOfRegion = EstimateModeComplex(coverage[i,mixOfSamples])
  modeOfRegion
}


for (i in 1:length(potentiallyPolymorphic)) {
  coverage.normalised[potentiallyPolymorphic[i]] = coverage.normalised[potentiallyPolymorphic[i]] / mediansOfPolymorphic[i]
}
