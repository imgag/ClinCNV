
if (!exists("polymorphicRegions")) {
  print("You need to specify a BED file with polymorphic regions you want to genotype. Genotyping without is not possible. QUit.")
  quit()
}

setwd(opt$folderWithScript)
source(paste0(opt$folderWithScript, "/germline/mCNVhelpers.R"))
setwd(opt$out)


if (!dir.exists(paste0(opt$out, "/normal/"))) {
  dir.create(paste0(opt$out, "/normal/"))
}
folder_name_mcnv <- paste0(opt$out, "/normal/mCNVs_genotyped/")

if (!dir.exists(folder_name_mcnv)) {
  dir.create(folder_name_mcnv)
}


vect_of_norm_likeliks = fast_dt_list(100)

EstimateModeSimple <- function(x) {
  density_of_x <-  density(x, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  mu
}


EstimateModeSimpleCov <- function(x) {
  tmpX = x[which(x > 0.25)]
  if (length(tmpX) < 10) {
    tmpX = x[which(x > quantile(x, 0.95))]
    if (median(tmpX) < 0.25) {
      return(1)
    }
    return(median(tmpX))
  }
  density_of_x <-  density(tmpX, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  mu
}


## MAKING COHORT ONLY OF LOW VARIANCE SAMPLES
QnSample <- apply(coverage[which(!bedFile[,1] %in% c("chrX","chrY")),], 2, mad)
modeOfVariances <- EstimateModeSimple(QnSample)
samplesForExclusion <- which(QnSample > 1.5 * modeOfVariances)
if (length(samplesForExclusion) > 0) {
  print(paste("Samples", paste(colnames(coverage)[samplesForExclusion], collapse=", "), "have variances 2.5 times bigger than the mode variance of the cohort. We have to exclude them."))
  coverageForNormalization = coverage[,-samplesForExclusion]
  genderOfSamplesUsed = genderOfSamples[-samplesForExclusion]
  if (ncol(coverageForNormalization) < 50) {
    paste("The amount of samples is too small for polymorphic calling!")
  }
} else {
  coverageForNormalization = coverage
  genderOfSamplesUsed = genderOfSamples
}


autosomes = which(!bedFile[,1] %in% c("chrX", "chrY"))
mediansAndSdsPolymorphic = calculateLocationAndScale(bedFile, coverageForNormalization, genderOfSamplesUsed, autosomes, polymorphic=T)
mediansOfPolymorphic = mediansAndSdsPolymorphic[[1]][,1]
coverage.normalised.polymorph = sweep(coverageForNormalization, 1, mediansOfPolymorphic + 10**-40, FUN="/")
rm(mediansAndSdsPolymorphic)
rm(coverageForNormalization)


regionsToRemove <- which(mediansOfPolymorphic <= 0.25 | mediansOfPolymorphic >= sqrt(16/2))
regionsToRemove = unique(c(regionsToRemove, regionsToRemove + 1, regionsToRemove - 1))
if (length(regionsToRemove) > 0) {
  bedFilePolymorph = bedFile[-regionsToRemove,]
  coverage.normalised.polymorph = coverage.normalised.polymorph[-regionsToRemove,]
  mediansOfPolymorphic = mediansOfPolymorphic[-regionsToRemove]
}

QnSample <- apply(coverage.normalised.polymorph[which(!bedFilePolymorph[,1] %in% c("chrX","chrY")),], 2, Sn)


matrixOfMultipliers <- matrix(0, nrow=1000, ncol=ncol(coverage.normalised.polymorph))
for (i in 1:nrow(matrixOfMultipliers)) {
  vec <- rnorm(ncol(coverage.normalised.polymorph), sd=QnSample)
  res = sd(vec)
  matrixOfMultipliers[i,] = QnSample / res
}

multipliersSamples <- apply(matrixOfMultipliers, 2, median)
sdsOfProbes <- apply(coverage.normalised.polymorph, 1, Qn)
bandwidths <- parApply(cl=cl, coverage.normalised.polymorph, 1, bw.SJ)
dataToTrain = data.frame(cbind(sdsOfProbes, bandwidths))
newlm <- rlm(sdsOfProbes ~ bandwidths, data=dataToTrain)
predictedVariances = predict(newlm, dataToTrain, interval="prediction", level=0.99)
sdsOfProbesCorrected = apply(cbind(sdsOfProbes, predictedVariances[,3]), 1, min)

# lower bound of SD
lowerBoundOfSD = quantile(sdsOfProbesCorrected, 0.0001)

copyNumberForReportingGlobal = matrix(0, nrow=nrow(polymorphicRegions), ncol(coverage.normalised.polymorph))
copyNumberForReportingGlobal = cbind(polymorphicRegions, copyNumberForReportingGlobal)
colnames(copyNumberForReportingGlobal) = c("chr", "start", "end", colnames(coverage.normalised.polymorph))


locations = list()
for (i in 1:16) {
  locations[[i]] = sqrt(1:20/i) #[sqrt(1:20/i) > 0.4]
}

if (!is.null(polymorphicRegions)) {
  for (region in 1:nrow(copyNumberForReportingGlobal)) {
    chrom  = polymorphicRegions[region,1]
  samplesForAnalysis = 1:ncol(coverage.normalised.polymorph)
    if (chrom == "chrY") next
     multipliersSamplesForAnalysis = multipliersSamples[samplesForAnalysis]
  
    
      whichInsideVariant <- which(bedFilePolymorph[,1]==chrom & as.numeric(bedFilePolymorph[,2]) >= as.numeric(polymorphicRegions[region,2]) - 10 & 
                                    as.numeric(bedFilePolymorph[,3]) <= as.numeric(polymorphicRegions[region,3]) + 10)
       if (length(whichInsideVariant) > 0) {
         mcnvCopyNumber <- findFinalState(coverage.normalised.polymorph[whichInsideVariant,samplesForAnalysis,drop=F], 
                                          bedFilePolymorph[whichInsideVariant,,drop=F],
                                          multipliersSamplesForAnalysis, cluster, T, F, folder_name_mcnv, 
                                          median(mediansOfPolymorphic[whichInsideVariant]))
          copyNumberForReportingGlobal[region,4:ncol(copyNumberForReportingGlobal)] = mcnvCopyNumber
       } else {
         print(region)
         copyNumberForReportingGlobal[region,4:ncol(copyNumberForReportingGlobal)] = "NA"
       }
  }
}

write.table(copyNumberForReportingGlobal, file=paste0(opt$out, "/normal/mCNVs_genotyped/", cluster, "mCNVs.txt"), quote = F, row.names = F, col.names = T, sep="\t")



