setwd(opt$folderWithScript)
source(paste0(opt$folderWithScript, "/germline/mCNVhelpers.R"))

load("/Users/gdemidov/Downloads/prepared-2.RData")



# signalToNoise = mediansAndSds[[2]][which(!bedFile[,1] %in% c("chrX", "chrY")),1] / (mediansAndSds[[2]][which(!bedFile[,1] %in% c("chrX", "chrY")),2])
# signalToNoiseSex = mediansAndSds[[2]][which(bedFile[,1] %in% c("chrX", "chrY")),1] / (mediansAndSds[[2]][which(bedFile[,1] %in% c("chrX", "chrY")),2])
# signalToNoiseAll = mediansAndSds[[2]][,1] / mediansAndSds[[2]][,2]
# potentiallyPolymorphic <- which(mediansAndSds[[2]][,1] > sqrt(4/2) | 
#                                   (signalToNoiseAll < quantile(signalToNoise, 0.05) & !(bedFile[,1] %in% c("chrX", "chrY"))) |
#                                   (signalToNoiseAll < quantile(signalToNoiseSex, 0.05) & bedFile[,1] %in% c("chrX", "chrY"))
#                                   )
# 
# potentiallyPolymorphic = unique(union(potentiallyPolymorphic, c(potentiallyPolymorphic + 1, potentiallyPolymorphic - 1)))
# potentiallyPolymorphic = potentiallyPolymorphic[which(potentiallyPolymorphic > 0 & potentiallyPolymorphic < nrow(coverage.normalised))]
# 
# potentiallyPolymorphicCoverageNormalised = coverage.normalised[potentiallyPolymorphic,]

  
mediansAndSdsPolymorphic = calculateLocationAndScale(bedFile, coverage, genderOfSamples, autosomes, T)
coverage.normalised.polymorph = mediansAndSdsPolymorphic[[1]]


mediansOfPolymorphic = mediansAndSdsPolymorphic[[2]][,1]
regionsToRemove <- which(mediansOfPolymorphic <= 0.1 | mediansOfPolymorphic >= sqrt(20/2))
if (length(regionsToRemove) > 0) {
  bedFilePolymorph = bedFile[-regionsToRemove,]
  coverage.normalised.polymorph = coverage.normalised.polymorph[-regionsToRemove,]
  mediansOfPolymorphic = mediansOfPolymorphic[-regionsToRemove]
}


