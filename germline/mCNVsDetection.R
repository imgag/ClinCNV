setwd(opt$folderWithScript)
source(paste0(opt$folderWithScript, "/germline/mCNVhelpers.R"))



vect_of_norm_likeliks = fast_dt_list(100)
  
autosomes = which(!bedFile[,1] %in% c("chrX", "chrY"))
mediansAndSdsPolymorphic = calculateLocationAndScale(bedFile, coverage, genderOfSamples, autosomes, T)
coverage.normalised.polymorph = mediansAndSdsPolymorphic[[1]]


mediansOfPolymorphic = mediansAndSdsPolymorphic[[2]][,1]
regionsToRemove <- which(mediansOfPolymorphic <= 0.25 | mediansOfPolymorphic >= sqrt(16/2))
regionsToRemove = unique(c(regionsToRemove, regionsToRemove + 1, regionsToRemove - 1))
if (length(regionsToRemove) > 0) {
  bedFilePolymorph = bedFile[-regionsToRemove,]
  coverage.normalised.polymorph = coverage.normalised.polymorph[-regionsToRemove,]
  mediansOfPolymorphic = mediansOfPolymorphic[-regionsToRemove]
}

QnSample <- apply(coverage.normalised.polymorph[which(!bedFilePolymorph[,1] %in% c("chrX","chrY")),], 2, Qn)
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
predictedVariances = predict(newlm, dataToTrain, interval="prediction", level=0.999)
sdsOfProbesCorrected = apply(cbind(sdsOfProbes, predictedVariances[,3]), 1, min)

# lower bound of SD
lowerBoundOfSD = quantile(sdsOfProbesCorrected, 0.0001)
vect_of_norm_likeliks <- fast_dnorm_list()
minimum_length_of_CNV = 1
threshold = 20
initial_state = 1
price_per_tile = 5
resultingMCNVs <- matrix(0, nrow=0, ncol=11)
copyNumberForReporting = matrix(0, nrow=0, ncol=3 + ncol(coverage))
colnames(copyNumberForReporting) = c("chr", "start", "end", colnames(coverage))
for (l in 1:length(left_borders)) {
  chrom = names(left_borders)[l]
  if (chrom == "chrX") next
  if (chrom == "chrY") next
  for (k in 1:2) {
    print(paste("Polymorphic calling at chrom", chrom, "and arm", k, "started"))
    which_to_allow <- "NA"
    which_to_allow_ontarget <- "NA"
    if (k == 1) {
      which_to_allow = which(bedFilePolymorph[,1] == chrom & as.numeric(bedFilePolymorph[,2]) <= as.numeric(left_borders[[l]]) )
    } else {
      which_to_allow = which(bedFilePolymorph[,1] == chrom & as.numeric(bedFilePolymorph[,2]) >= as.numeric(right_borders[[l]]) )
    }
    if (length(which_to_allow) <= 1) {
      next
    }
    
    coverageToWorkWith = coverage.normalised.polymorph[which_to_allow,]
    mediansOfPolymorphicLocal = mediansOfPolymorphic[which_to_allow]
    localSdsOfProbes <- sdsOfProbesCorrected[which_to_allow]
    toyBedFilePolymorph = bedFilePolymorph[which_to_allow,]
    locations = list()
    for (i in 1:16) {
      locations[[i]] = sqrt(1:20/i)[sqrt(1:20/i) > 0.4]
    }
    matrixOfLikeliksLocal = matrix(threshold, nrow= nrow(coverageToWorkWith), ncol=2)
    matrixOfLikeliksLocal[,1] = 0
    correlations <- sapply(2:nrow(coverageToWorkWith), function(i) {
      moreThanZero = which(coverageToWorkWith[i,] > 0.25 & coverageToWorkWith[i-1,] > 0.25) ; 
      cor(coverageToWorkWith[i,moreThanZero], coverageToWorkWith[i-1,moreThanZero])})
    for (i in 1:nrow(coverageToWorkWith)) {
      correlationsAround = c(correlations[i - 1])
      if (i < nrow(coverageToWorkWith) - 1) {
        correlationsAround = c(correlationsAround, correlations[i])
      }
      if (i %% 1000 == 0) print(i)
      coverageOfProbe = coverageToWorkWith[i,]
      if (max(correlationsAround) > quantile(correlations, 0.9)) {
        notHomozygousDeletions = which(coverageOfProbe > 0.25)
        if (length(notHomozygousDeletions) < length(coverageOfProbe)) {
          homozygousDelShit = median(coverageOfProbe[which(coverageOfProbe < 0.25)]) ** 2
          coverageOfProbe = sqrt(abs(coverageOfProbe ** 2 - homozygousDelShit))
        }
        sdNormalised = localSdsOfProbes[i]
        matrixOfLikeliksLocal[i,1] = -2 * likelihoodOfGaussianMixture(c(1), coverageOfProbe[notHomozygousDeletions], sdNormalised , 
                                                                 0.05 * (length(notHomozygousDeletions)) / ncol(coverageToWorkWith), 0.1, 
                                                                 multipliersSamples[notHomozygousDeletions], lowerBoundOfSD)
        for (j in 1:length(locations)) {
          if (mediansOfPolymorphicLocal[i] < sqrt(1/2) & j >= 3) next
          if (mediansOfPolymorphicLocal[i] < sqrt(3/2) & j >= 6) next
          if (mediansOfPolymorphicLocal[i] > sqrt(4/2) & j <= 2) next
          if (mediansOfPolymorphicLocal[i] > sqrt(8/2) & j <= 4) next
          sdNormalised = localSdsOfProbes[i]
          vecOfMeans = locations[[j]]
          vecOfMeans = vecOfMeans[which(vecOfMeans > min(coverageOfProbe[notHomozygousDeletions]) & vecOfMeans < max(coverageOfProbe[notHomozygousDeletions]))]
          if (length(vecOfMeans) == 1) next
          likelikAndWeights = likelihoodOfGaussianMixture(vecOfMeans, coverageOfProbe[notHomozygousDeletions], sdNormalised, 
                                                          0.05 * (length(notHomozygousDeletions)) / ncol(coverageToWorkWith), 0.1, 
                                                          multipliersSamples[notHomozygousDeletions], lowerBoundOfSD)
          tmpLikelik =  (
            -2 * likelikAndWeights[[1]] + (length(which(likelikAndWeights[[2]] > 0.1 / length(notHomozygousDeletions)))  + 1) * log(length(notHomozygousDeletions))
                           )
          if (tmpLikelik < matrixOfLikeliksLocal[i,2] & tmpLikelik < matrixOfLikeliksLocal[i,1]) {
            matrixOfLikeliksLocal[i,2] = tmpLikelik
            localSdsOfProbes[i] = likelikAndWeights[[3]]
          }
        }
      } else {
        matrixOfLikeliksLocal[i,2] = threshold
      }
    }
    
    found_CNVsTmp <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state,  matrixOfLikeliksLocal, initial_state))

    finalMCNVs <- matrix(0, nrow=0, ncol=ncol(found_CNVsTmp))
    smallRegions <- matrix(0, nrow=0, ncol=ncol(found_CNVsTmp))
    regionsToLookForMCNVs <- found_CNVsTmp
    if (nrow(found_CNVsTmp) > 0) {
      counter = 1
      while (counter <= nrow(regionsToLookForMCNVs)) {
        if (regionsToLookForMCNVs[counter,3] - regionsToLookForMCNVs[counter,2] < 2) {
          smallRegions = rbind(smallRegions, regionsToLookForMCNVs[counter,])
        } else if (regionsToLookForMCNVs[counter,3] - regionsToLookForMCNVs[counter,2] < 4) {
          finalMCNVs = rbind(finalMCNVs, regionsToLookForMCNVs[counter,])
        } else {
          wasDivided = F
          for (j in (regionsToLookForMCNVs[counter,2] + 1):(regionsToLookForMCNVs[counter,3] - 2)) {
            if (!checkConnectivity(coverageToWorkWith[j,], coverageToWorkWith[j + 1,])) {
              if (j + 2 <= nrow(coverageToWorkWith)) {
                if (checkConnectivity(coverageToWorkWith[j,], coverageToWorkWith[j + 2,])) {
                  next
                }
              }
              leftVar = regionsToLookForMCNVs[counter,]
              leftVar[3] = j
              rightVar = regionsToLookForMCNVs[counter,]
              rightVar[2] = j + 1
              regionsToLookForMCNVs = rbind(regionsToLookForMCNVs, leftVar)
              regionsToLookForMCNVs = rbind(regionsToLookForMCNVs, rightVar)
              wasDivided = T
              break
            }
          }
          if (!wasDivided) {
            finalMCNVs = rbind(finalMCNVs, regionsToLookForMCNVs[counter,])
          }
        }
        counter = counter + 1
      }
    }
    finalMCNVs = rbind(finalMCNVs, smallRegions)
    
    if (nrow(finalMCNVs) > 0)
    for (i in 1:nrow(finalMCNVs)) {
      mcnvCopyNumber <- findFinalState(coverageToWorkWith[finalMCNVs[i,2]:finalMCNVs[i,3],], 
                     median(mediansOfPolymorphicLocal[(finalMCNVs[i,2] + 1):(finalMCNVs[i,3] - 1)]), 
                     toyBedFilePolymorph[finalMCNVs[i,2]:finalMCNVs[i,3],],
                     multipliersSamples)
      if (length(which(mcnvCopyNumber != 2)) < 0.05 * ncol(coverageToWorkWith)) {
        print(i)
        next
      }
      coverageNeededToCheck = coverageToWorkWith[finalMCNVs[i,2]:finalMCNVs[i,3],]
      medianOfCoverage = median(mediansOfPolymorphicLocal[(finalMCNVs[i,2] + 1):(finalMCNVs[i,3] - 1)])
      sdNormalised = median(localSdsOfProbes[(finalMCNVs[i,2] + 1):(finalMCNVs[i,3] - 1)])
      toyBedFilePolymorphCurrent = toyBedFilePolymorph[finalMCNVs[i,2]:finalMCNVs[i,3],]
      copyNumberForReporting = rbind(copyNumberForReporting, c(toyBedFilePolymorph[finalMCNVs[i,2],1], 
                                                               toyBedFilePolymorph[finalMCNVs[i,2],2],
                                                               toyBedFilePolymorph[finalMCNVs[i,3],3],
                                                               mcnvCopyNumber))
    }
   }
}

folder_name <- paste0(opt$out, "/normal/")
write.table(copyNumberForReporting, file=paste0(opt$out, "/normal/", cluster, "mCNVs.txt"), quote = F, row.names = F)
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}
for (sam_no in 4:ncol(copyNumberForReporting)) {
  sample_name <- colnames(copyNumberForReporting)[sam_no]
  print(paste("Working with germline sample", sample_name, Sys.time()))
  if (!dir.exists(paste0(folder_name, sample_name))) {
    dir.create(paste0(folder_name, sample_name))
  }
  makeTrackAnnotation <- function(fileName) {
    if (!file.exists(fileName)) {
      file.create(fileName)
      fileConn<-file(fileName)
      writeLines(c("#type=GENE_EXPRESSION",
                   paste0("#track graphtype=points name=\"", paste(sample_name, "polymorphic"), "\" color=0,0,255 altColor=255,0,0 maxHeightPixels=80:80:80 viewLimits=0:2:6 yLineMark=2 yLineOnOff=on"),
                   paste("ID", "chr", "start", "end", "CN", "loglik", "value", sep="\t")), fileConn)
      
      close(fileConn)
    }
  }
  fileToOut <- paste0(folder_name, sample_name, paste0("/", sample_name, "_mcnvs.seg"))
  makeTrackAnnotation(fileToOut)
  for (i in 1:nrow(copyNumberForReporting)) {
    detectedCN = as.numeric(copyNumberForReporting[i,sam_no])
    copyNumberSegment = matrix(c(sample_name, copyNumberForReporting[i,1], copyNumberForReporting[i,2], copyNumberForReporting[i,3], detectedCN, 1000, detectedCN), nrow=1, ncol=7)
    write(paste(copyNumberSegment[1,], collapse="\t"), file=fileToOut, append=TRUE)
  }
}
