setwd(opt$folderWithScript)
source(paste0(opt$folderWithScript, "/germline/mCNVhelpers.R"))

load("/Users/gdemidov/Downloads/prepared-2.RData")



signalToNoise = mediansAndSds[[2]][which(!bedFile[,1] %in% c("chrX", "chrY")),1] / (mediansAndSds[[2]][which(!bedFile[,1] %in% c("chrX", "chrY")),2])
signalToNoiseSex = mediansAndSds[[2]][which(bedFile[,1] %in% c("chrX", "chrY")),1] / (mediansAndSds[[2]][which(bedFile[,1] %in% c("chrX", "chrY")),2])
signalToNoiseAll = mediansAndSds[[2]][,1] / mediansAndSds[[2]][,2]
potentiallyPolymorphic <- which(mediansAndSds[[2]][,1] > sqrt(4/2) | 
                                  (signalToNoiseAll < quantile(signalToNoise, 0.05) & !(bedFile[,1] %in% c("chrX", "chrY"))) |
                                  (signalToNoiseAll < quantile(signalToNoiseSex, 0.05) & bedFile[,1] %in% c("chrX", "chrY"))
                                  )

potentiallyPolymorphic = unique(union(potentiallyPolymorphic, c(potentiallyPolymorphic + 1, potentiallyPolymorphic - 1)))
potentiallyPolymorphic = potentiallyPolymorphic[which(potentiallyPolymorphic > 0 & potentiallyPolymorphic < nrow(coverage.normalised))]

potentiallyPolymorphicCoverageNormalised = coverage.normalised[potentiallyPolymorphic,]
mediansOfPolymorphic = foreach (i=1:length(potentiallyPolymorphic), .combine="c") %dopar% {
  modeOfRegion = EstimateModeSimple(potentiallyPolymorphicCoverageNormalised[i,])
  modeOfRegion
}
rm(potentiallyPolymorphicCoverageNormalised)

regionsToRemove <- potentiallyPolymorphic[which(mediansOfPolymorphic == 0)]
if (length(regionsToRemove) > 0) {

bedFile = bedFile[-regionsToRemove,]
coverage.normalised = coverage.normalised[-regionsToRemove,]
potentiallyPolymorphic = potentiallyPolymorphic[-which(mediansOfPolymorphic == 0)]
}

for (i in 1:length(potentiallyPolymorphic)) {
  coverage.normalised[potentiallyPolymorphic[i],] = coverage.normalised[potentiallyPolymorphic[i],] / mediansOfPolymorphic[i]
}



