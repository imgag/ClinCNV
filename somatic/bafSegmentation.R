
source(paste0(opt$folderWithScript, "/somatic/helpersBalleleFreq.R"))
setwd(opt$bafFolder)
clusterExport(cl, c('determineHeterozygousPositions', 'determineHeterozygousPositionsOverdispersed'))

makeTrackAnnotation <- function(fileName, ID, viewLimits, trackType="points", addText="", color="0,0,255") {
  file.create(fileName)
  fileConn<-file(fileName)
  writeLines(c("#type=GENE_EXPRESSION",
               paste0("#track graphtype=", trackType, " name=\"", ID, "\" color=", color, " altColor=255,0,0 maxHeightPixels=80:80:80 viewLimits=", viewLimits, " ", addText)), fileConn)
  close(fileConn)
}


returnBAlleleFreqs <- function(healthySampleName, tumorSampleName, folderBAF, bedFileForFiltering) {
  setwd(folderBAF)
  healthySample = NULL
  tumorSample = NULL
  if (file.exists(paste0(healthySampleName, ".tsv")) & file.exists(paste0(tumorSampleName, ".tsv"))) {
    # CHECK IF FILE IS EMPTY
    info = file.info(paste0(healthySampleName, ".tsv"), paste0(tumorSampleName, ".tsv"))
    empty = rownames(info[info$size == 0, ])
    
    if (length(empty) > 0) {
      return(list(NULL, NULL))
    }
    
    # FILE IS NOT EMPTY => READING
    healthySample <- read.table(paste0(healthySampleName, ".tsv"), stringsAsFactors = F, header=F, sep="\t")
    if (!startsWith(healthySample[1,1], "chr")) {
      healthySample[,1] = paste0("chr", healthySample[,1])
    }
    
    tumorSample <- read.table(paste0(tumorSampleName, ".tsv"), stringsAsFactors=F, header=F, sep="\t")
    if (!startsWith(tumorSample[1,1], "chr")) {
      tumorSample[,1] = paste0("chr", tumorSample[,1])
    }
    if (nrow(healthySample) == 0) {
      return(list(NULL, NULL))
    }
    whichAreNAHealthy = which(healthySample[,5] == "n/a")
    whichAreNATumor = which(tumorSample[,5] == "n/a")
    
    if (length(whichAreNAHealthy) > 0) {
      healthySample = healthySample[-whichAreNAHealthy,]
    }
    if (length(whichAreNATumor) > 0) {
      tumorSample = tumorSample[-whichAreNATumor,]
    }
    
    if (ncol(healthySample) == 5) {
      healthySample = cbind(healthySample[,1:3], apply(healthySample[,1:3], 1, function(x) {paste0(x, collapse="_")}), healthySample[,4:5])
    }
    if (ncol(tumorSample) == 5) {
      tumorSample = cbind(tumorSample[,1:3], apply(tumorSample[,1:3], 1, function(x) {paste0(x, collapse="_")}), tumorSample[,4:5])
    }
    
    colnames(healthySample) <- c("chr", "start", "end", "Feature", "freq", "depth")
    colnames(tumorSample) <- c("chr", "start", "end", "Feature", "freq", "depth")
    
    healthySample = healthySample[which(healthySample[,6] > max(median(healthySample[,6]) / 10, 30) & healthySample[,6] < quantile(healthySample[,6], 0.975)),]
    tumorSample = tumorSample[which(tumorSample[,6] > max(median(tumorSample[,6]) / 10, 30)),]
    
    indicesOfSNVsToRemove <- c()
    currentChrom = "chrN"
    matrOfBedRegionsInChrom = c()
    for (i in 1:nrow(healthySample)) {
      if (healthySample[i,1] != currentChrom) {
        currentChrom = healthySample[i,1]
        matrOfBedRegionsInChrom = bedFileForFiltering[which(bedFileForFiltering[,1] == currentChrom),]
      }
      ifItIsInsideBed <- which(matrOfBedRegionsInChrom[,2] - 10 <= healthySample[i,2] & matrOfBedRegionsInChrom[,3] + 10 >= healthySample[i,3])
      if (length(ifItIsInsideBed) == 0) {
        indicesOfSNVsToRemove <- c(indicesOfSNVsToRemove, i)
      }
    }
    if (length(indicesOfSNVsToRemove) > 0)
      healthySample = healthySample[-indicesOfSNVsToRemove,]
    
    i = 1
    indicesToRemove <- c()
    lengthOfClusteredVariants <- 30
    while (i != nrow(healthySample) - 1) {
      if (abs(healthySample[i,2] - healthySample[i + 1,2]) < lengthOfClusteredVariants) {
        indicesToRemove <- c(indicesToRemove, i)
        indicesToRemove <- c(indicesToRemove, i + 1)
      }
      i = i + 1
    }
    if (length(indicesToRemove) > 0)
      healthySample = healthySample[-unique(indicesToRemove),]
    
    
    
    healthyFeatures <- healthySample$Feature
    tumorFeatures <- tumorSample$Feature
    healthyFeatures = healthyFeatures[!(duplicated(healthyFeatures) | duplicated(healthyFeatures, fromLast = TRUE)) ]
    tumorFeatures = tumorFeatures[!(duplicated(tumorFeatures) | duplicated(tumorFeatures, fromLast = TRUE)) ]
    featuresPresentedInBoth <- intersect(healthyFeatures, tumorFeatures)
    
    tumorSample = tumorSample[which(tumorSample$Feature %in% featuresPresentedInBoth),]
    healthySample = healthySample[which(healthySample$Feature %in% featuresPresentedInBoth),]
    
    tumorSample = tumorSample[order(tumorSample[,1],tumorSample[,2] ),]
    healthySample = healthySample[order(healthySample[,1], healthySample[,2]),]
    
    # Determining positions which are heterozygous 
    potentiallyHeterozygous = as.numeric(healthySample[which(!(healthySample[,1] %in% c("chrX","chrY","X","Y")) & as.numeric(healthySample[,6]) > 30),5])
    potentiallyHeterozygous <- potentiallyHeterozygous[which(potentiallyHeterozygous > 0.25 & potentiallyHeterozygous < 0.75)]
    print(summary(potentiallyHeterozygous))
    if (!is.na(potentiallyHeterozygous)) {
      heterozygousAlleleShift <- median(potentiallyHeterozygous)
    } else { heterozygousAlleleShift = 0.48 }
    print(paste("Heterozygous allele shift for particular sample", heterozygousAlleleShift))
    overdispersionCorrectionNormal = extractVariancesFromBAF(healthySample, heterozygousAlleleShift)
    #overdispersionCorrectionTumor = extractVariancesFromBAF(tumorSample, heterozygousAlleleShift)
    if (is.null(overdispersionCorrectionNormal)) {
      heterozygousPositions <- parApply(cl=cl, healthySample[,5:6], 1, function(vec) {determineHeterozygousPositions(as.numeric(vec[1]), as.numeric(vec[2]), heterozygousAlleleShift)})
    } else {
      overdispersionFactors = rep(0, nrow(healthySample[,5:6]))
      for (l in 1:nrow(healthySample)) {
        closestDepth = which.min(abs(healthySample[l,6] - overdispersionCorrectionNormal[,2]))
        overdispersionFactors[l] = overdispersionCorrectionNormal[closestDepth,1]
      }
      heterozygousPositions <- parApply(cl=cl, healthySample[,5:6], 1, function(vec) {determineHeterozygousPositionsOverdispersed(as.numeric(vec[1]), as.numeric(vec[2]), heterozygousAlleleShift, overdispersionFactors)})
    }
    healthySample = healthySample[heterozygousPositions, ]
    tumorSample = tumorSample[heterozygousPositions, ]
    
    meaningfulResult <- list()
    meaningfulResult[[healthySampleName]] = healthySample
    meaningfulResult[[tumorSampleName]] = tumorSample
    return(meaningfulResult)
  } else {
    return(list(healthySample, tumorSample))
  }
}


determineAllowedChroms <- function(healthySample, tumorSample, healthySampleName, tumorSampleName, folderBAF, leftBorders, rightBorders, endsOfChroms) {
  # OVERDISPERDSION CORRECTION
  heterozygousAlleleShift = median(as.numeric(healthySample[,5]))
  overdispersionCorrectionNormal = extractVariancesFromBAF(healthySample, heterozygousAlleleShift)
  overdispersionCorrectionTumor = extractVariancesFromBAF(tumorSample, heterozygousAlleleShift)
  vecOfPvalues <- rep(0, nrow(tumorSample))
  if (!is.null(overdispersionCorrectionNormal) & !is.null(overdispersionCorrectionTumor)) {
    overdispersionFactorsNornm = rep(0, nrow(healthySample[,5:6]))
    overdispersionFactorsTum = rep(0, nrow(tumorSample[,5:6]))
    for (l in 1:nrow(healthySample)) {
      closestDepthNorm = which.min(abs(healthySample[l,6] - overdispersionCorrectionNormal[,2]))
      overdispersionFactorsNornm[l] = overdispersionCorrectionNormal[closestDepthNorm,1]
      closestDepthTum = which.min(abs(tumorSample[l,6] - overdispersionCorrectionTumor[,2]))
      overdispersionFactorsTum[l] = overdispersionCorrectionTumor[closestDepthTum,1]
      
      refAleleTum <- (round(as.numeric(tumorSample[l,5]) * as.numeric(tumorSample[l,6])))
      altAleleTum <- as.numeric(tumorSample[l,6]) - refAleleTum
      refAleleNorm<- (round(as.numeric(healthySample[l,5]) * as.numeric(healthySample[l,6])))
      altAleleNorm <- as.numeric(healthySample[l,6]) - refAleleNorm
      
      vecOfPvalues[l] = passPropTestVarCorrection(refAleleTum, refAleleNorm, altAleleTum, altAleleNorm, overdispersionFactorsNornm[l], overdispersionFactorsTum[l])
    }
  } else {
    for (i in 1:nrow(tumorSample)) {
      refAleleTum <- (round(as.numeric(tumorSample[i,5]) * as.numeric(tumorSample[i,6])))
      altAleleTum <- as.numeric(tumorSample[i,6]) - refAleleTum
      refAleleNorm<- (round(as.numeric(healthySample[i,5]) * as.numeric(healthySample[i,6])))
      altAleleNorm <- as.numeric(healthySample[i,6]) - refAleleNorm
      vecOfPvalues[i] = passPropTest(refAleleTum, refAleleNorm, altAleleTum, altAleleNorm)
    }
  }
  
  pvalues <- round(vecOfPvalues, digits=4)
  pvalues <- format(pvalues, scientific = FALSE)
  pvalues <- matrix(pvalues, ncol=1)
  values <- rep(1, length(pvalues))
  values[which(pvalues < 0.05)] = 0
  values[which(pvalues < 0.01)] = -1
  values = matrix(values, ncol=1)
  trackToWriteOut <- cbind(healthySample[,c(1,2,3)], pvalues, values)
  colnames(trackToWriteOut)[4] = "Features"
  
  thresholdOfNonNormalVariant = 0.01 # pvalue, if lower = BAF is different in normal and tumor
  chroms = paste0("chr", 1:22)
  chroms = c(chroms, "chrX")
  numberOfSNVs = rep(0, 2 * length(left_borders))
  evaluated = rep(0, 2 * length(left_borders))
  namesOfChromArms <- rep(0, 2 * length(left_borders))
  plotLabels <- rep(0, 2 * length(left_borders))
  counter = 1
  for (l in 1:length(left_borders)) {
    chrom = names(left_borders)[l]
    start = left_borders[[l]]
    end = right_borders[[l]]
    endOfChrom = endsOfChroms[[l]]
    for (k in 1:2) {
      if (k == 1) {
        rowsFromChrom = which(healthySample[,1] == chrom & healthySample[,2] <= as.numeric(start) )
        namesOfChromArms[counter] = paste(chrom, 0, start, sep="-")
        plotLabels[counter] = paste(chrom, "left #", length(rowsFromChrom), "SNVs")
      } else {
        rowsFromChrom = which(healthySample[,1] == chrom & healthySample[,2] >= as.numeric(end) )
        namesOfChromArms[counter] = paste(chrom, end, endOfChrom, sep="-")
        plotLabels[counter] = paste(chrom, "right #", length(rowsFromChrom), "SNVs")
      }
      
      if (length(rowsFromChrom) > 15) {
        evaluationOfChorm = (length(which(pvalues[rowsFromChrom] < thresholdOfNonNormalVariant))) / (length(rowsFromChrom))
      } else {
        evaluationOfChorm = 1.0
      }
      
      evaluated[counter] = evaluationOfChorm
      numberOfSNVs[counter] = length(rowsFromChrom)
      
      counter = counter + 1
    }
  }
  indicesOfAllowedChroms = which(evaluated < max(sort(evaluated)[3], 0.05))
  colVec <- rep("red", length(evaluated))
  indicesOfAllowedButNotBestChroms = which(evaluated > 0.05 & evaluated < sort(evaluated)[6])
  colVec[indicesOfAllowedChroms] = "darkgreen"
  colVec[indicesOfAllowedButNotBestChroms] = "darkorange"
  names(evaluated) = namesOfChromArms
  allowedChroms = names(evaluated)[union(indicesOfAllowedChroms, indicesOfAllowedButNotBestChroms)]
  
  
  
  subDir = paste0(tumorSampleName, "_", healthySampleName)
  dir.create(file.path(folderBAF, "result", subDir))
  setwd(file.path(folderBAF, "result", subDir))
  png(paste0(tumorSampleName, "_", healthySampleName, ".png"), width=1400, height=600)
  op <- par(mar=c(11,4,4,2))
  x <- barplot(evaluated, col=colVec, main=paste(tumorSampleName, healthySampleName), ylim=c(0,1), xaxt="n")
  text(x-2.5, par("usr")[3] - 0.15, labels = plotLabels, srt = 45, pos = 1, xpd = TRUE)
  par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  dev.off()
  
  png("overdispersion.png", width=1024, height=1024)
  plot(overdispersionFactorsNornm ~ overdispersionFactorsTum, main="Dispersion over Binomial")
  dev.off()
  
  trackFilename = paste0(tumorSampleName, "_", healthySampleName, ".igv")
  print(getwd())
  print(trackFilename)
  makeTrackAnnotation(trackFilename, trackFilename, "-1:1", "line")
  write.table(trackToWriteOut, file=trackFilename, quote = F, row.names = F, sep="\t", append=T)
  
  healthyBafFilename = paste0(healthySampleName, ".igv")
  makeTrackAnnotation(healthyBafFilename, healthySampleName, "0:1", "points", "yLineMark=1 yLineOnOff=on", "0,120,0")
  trackToWriteOut <- cbind(healthySample[,c(1,2,3,4,5)])
  write.table(trackToWriteOut, file=healthyBafFilename, quote = F, row.names = F, sep="\t", append=T)
  
  tumorBafFilename = paste0(tumorSampleName, ".igv")
  makeTrackAnnotation(tumorBafFilename, tumorSampleName, "0:1", "points", "yLineMark=1 yLineOnOff=on", "120,0,0")
  trackToWriteOut <- cbind(tumorSample[,c(1,2,3,4,5)])
  write.table(trackToWriteOut, file=tumorBafFilename, quote = F, row.names = F, sep="\t", append=T)
  
  return(list(allowedChroms, overdispersionFactorsNornm, overdispersionFactorsTum))
}


returnAllowedChromsBaf <- function(pairs, normalCov, tumorCov, inputFolderBAF, bedFileForFiltering, leftBorders, rightBorders, endsOfChroms) {
  allowedChromsBaf <- list()
  bAlleleFreqsAllSamples <- list()
  overdispersionNormal = list()
  overdispersionTumor = list()
  for (i in 1:ncol(normalCov)) {
    sampleNames1 <- which(pairs[,2] == colnames(normalCov)[i])
    for (sampleName1 in sampleNames1) {
      sampleName2 <- which(colnames(tumorCov) == pairs[sampleName1,1])
      bAlleleFreqs <- returnBAlleleFreqs(colnames(normalCov)[i], colnames(tumorCov)[sampleName2], inputFolderBAF, bedFileForFiltering)
      if (!is.null(bAlleleFreqs[[1]])) {
        listOfValues =determineAllowedChroms(bAlleleFreqs[[1]], bAlleleFreqs[[2]],
                                            colnames(normalCov)[i], colnames(tumorCov)[sampleName2],
                                            inputFolderBAF,
                                            leftBorders, rightBorders, endsOfChroms)
        allowedChromsBaf[[paste(colnames(tumorCov)[sampleName2], colnames(normalCov)[i], sep="-")]] = listOfValues[[1]]
        overdispersionNormal[[paste(colnames(tumorCov)[sampleName2], colnames(normalCov)[i], sep="-")]] = listOfValues[[2]]
        overdispersionTumor[[paste(colnames(tumorCov)[sampleName2], colnames(normalCov)[i], sep="-")]] = listOfValues[[3]]
        bAlleleFreqsAllSamples[[paste(colnames(tumorCov)[sampleName2], colnames(normalCov)[i], sep="-")]] = bAlleleFreqs
      } else {
        print("BAF did not work well")
        print(i)
        print(colnames(tumorCov)[sampleName2])
        print(colnames(normalCov)[i])
      }
    }
  }
  return(list(allowedChromsBaf, bAlleleFreqsAllSamples, overdispersionNormal, overdispersionTumor))
}

mergeClustersCloseToEachOther <- function(means, volumes, radius = 0.05) {
  smallVolumes <- which(volumes < 0.01)
  if (length(smallVolumes) > 0) {
    means = means[-smallVolumes]
    volumes = volumes[-smallVolumes]
    volumes = volumes / sum(volumes)
  }
  
  alreadyTaken = rep(F, length(means))
  newMeans <- c()
  newVolumes <- c()
  for (i in 1:length(means)) {
    if (alreadyTaken[i] == F) {
      whichAreAroundMeanI = which(abs(means - means[i]) < radius & alreadyTaken == F)
      alreadyTaken[whichAreAroundMeanI] = T
      newMeans <- c(newMeans, sum(volumes[whichAreAroundMeanI] * means[whichAreAroundMeanI] / sum(volumes[whichAreAroundMeanI])))
      newVolumes <- c(newVolumes, sum(volumes[whichAreAroundMeanI]))
    }
  }
  return(matrix(rbind(newMeans, newVolumes), nrow=2))
}


predictWholeGenomeEvent <- function(healthySampleBAF, tumorSampleBAF, matrixOfCoverages, allowedChromsBafForThatSamples, bedFile) {
  purDel = 1.0
  purDup = 1.0
  
  logFoldChanges <- log2(matrixOfCoverages[1,]) - log2(matrixOfCoverages[2,])
  
  shift <- median(logFoldChanges[which(bedFile[matrixOfCoverages[3,],1] %in% allowedChromsBafForThatSamples)])
  print(paste("SHIFT", shift))
  logFoldChanges = logFoldChanges - shift
  
  toFilter <- which(logFoldChanges < -1.0 | logFoldChanges == 0)
  
  values <- as.numeric(healthySampleBAF[,5]) - as.numeric(tumorSampleBAF[,5])
  
  if (length(toFilter) > 0) {
    logFoldChanges = logFoldChanges[-toFilter]
    values = values[-toFilter]
  }
  
  
  moreThan0 = which(logFoldChanges > 0)
  mod <- densityMclust(values[moreThan0], modelNames=c("E"))
  means = mod$parameters$mean
  means = means[which(means > 0.05 | means < -0.05)]
  if (length(means) > 0) {
    minAbs <- min(abs(means)) / 2
  } else {
    minAbs = 0.05
  }
  positiveValues <- values[moreThan0]
  positiveValues[which(positiveValues < -minAbs)] = abs(positiveValues[which(positiveValues < -minAbs)] )
  mod <- densityMclust(positiveValues, modelNames=c("E"))
  plot(mod, what="density", data=positiveValues, breaks=50,xlim=c(-minAbs,0.5))
  mergedData <- (mergeClustersCloseToEachOther(mod$parameters$mean, mod$parameters$pro))
  if (ncol(mergedData) > 1) {
    mergedData <- mergedData[,order(mergedData[1,])]
    if (sum(mergedData[2,2:col(mergedData)]) > 0.1) {
      print("For deletions:")
      print(mergedData)
      print("PURITY:")
      purDel = (estimatePurity("del", min(max(mergedData[1,]), 0.5)))
      if (purDel < 1.0)
        print(purDel)
      else 
        print(estimatePurity("del", mergedData[1,1]))
    }
    if (sum(mergedData[2,2:col(mergedData)]) > 0.2 & purDel > 0.5) {
      purDup = purDel
    }
  }
  #readline(prompt="Press [enter] to continue")
  
  lessThan0 = which(logFoldChanges < 0)
  mod <- densityMclust(values[lessThan0], modelNames=c("E"))
  means = mod$parameters$mean
  means = means[which(means > 0.05 | means < -0.05)]
  if (length(means) > 0) {
    minAbs <- min(abs(means)) / 2
  } else {
    minAbs = 0.05
  }
  positiveValues <- values[lessThan0]
  positiveValues[which(positiveValues < -minAbs)] = abs(positiveValues[which(positiveValues < -minAbs)] )
  mod <- densityMclust(positiveValues, modelNames=c("E"))
  plot(mod, what="density", data=positiveValues, breaks=50,xlim=c(-minAbs,0.5))
  mergedData <- (mergeClustersCloseToEachOther(mod$parameters$mean, mod$parameters$pro))
  if (ncol(mergedData) > 1) {
    mergedData <- mergedData[,order(mergedData[1,])]
    if (sum(mergedData[2,2:col(mergedData)]) > 0.1) {
      print("For duplications:")
      print(mergedData)
      print("PURITY:")
      purDup = (estimatePurity("dup", mergedData[1,2]))
      if (purDup < 1.0)
        print(purDup)
      else 
        print(estimatePurity("dup", mergedData[1,1]))
    }
  }
  
  
  #readline(prompt="Press [enter] to continue")
  
  print(paste("Max purity:", max(purDel, purDup)))
  
  return(min(1.0, max(purDel, purDup)))
}

extractCoveragesHavingListOfSNVs <- function(sampleBAF, sampleNormalCoverage, sampleTumorCoverage, bedFile) {
  coverageOfSampleWithinSNVRegionsNormal <- rep(1, nrow(sampleBAF))
  coverageOfSampleWithinSNVRegionsTumor <- rep(1, nrow(sampleBAF))
  modifiedBedFile <- rep(1, nrow(sampleBAF))
  for (i in 1:nrow(sampleBAF)) {
    coordInBed <- which(bedFile[,1] == sampleBAF[i,1] & bedFile[,2] <= sampleBAF[i,2] & bedFile[,3] >= sampleBAF[i,3])
    if (length(coordInBed) == 1) {
      coverageOfSampleWithinSNVRegionsNormal[i] = sampleNormalCoverage[coordInBed]
      coverageOfSampleWithinSNVRegionsTumor[i] = sampleTumorCoverage[coordInBed]
      modifiedBedFile[i] = coordInBed
    } 
    if (length(coordInBed) == 0) {
      bedFromSameChrom <- which(bedFile[,1] == sampleBAF[i,1])
      distancesLeft = abs(bedFile[bedFromSameChrom,2] - sampleBAF[i,2])
      distancesRight = abs(bedFile[bedFromSameChrom,3] - sampleBAF[i,3])
      minLeft = min(distancesLeft)
      minRight = min(distancesRight)
      if (minLeft < 10**4 | minRight < 10**4) {
        if (minLeft < minRight) {
          coordInBed = bedFromSameChrom[which.min(distancesLeft)]
        } else {
          coordInBed = bedFromSameChrom[which.min(distancesRight)]
        }
        coverageOfSampleWithinSNVRegionsNormal[i] = sampleNormalCoverage[coordInBed]
        coverageOfSampleWithinSNVRegionsTumor[i] = sampleTumorCoverage[coordInBed]
        modifiedBedFile[i] = coordInBed
      }
    }
  }
  return(rbind(coverageOfSampleWithinSNVRegionsNormal, coverageOfSampleWithinSNVRegionsTumor, modifiedBedFile))
}



returnPurityPloidy <- function(pairs, normalCov, tumorCov, inputFolderBAF, bedFile, allowedChromsBaf) {
  
  purityPloidy <- list()
  for (i in 1:ncol(normalCov)) {
    sampleNames1 <- which(pairs[,2] == colnames(normalCov)[i])
    for (sampleName1 in sampleNames1) {
      sampleName2 <- which(colnames(tumorCov) == pairs[sampleName1,1])
      bAlleleFreqs <- returnBAlleleFreqs(colnames(normalCov)[i], colnames(tumorCov)[sampleName2], inputFolderBAF, bedFile)
      if (!is.null(bAlleleFreqs[[1]])) {
        
        allowedChromsBafForThatSamples <- allowedChromsBaf[[paste(colnames(tumorCov)[sampleName2], colnames(normalCov)[i], sep="-")]]
        matrixOfCoverages <- extractCoveragesHavingListOfSNVs(bAlleleFreqs[[1]], normalCov[,i], tumorCov[,sampleName2], bedFile)
        if (length(bAlleleFreqs) >= 2 & !is.null(bAlleleFreqs[[1]])) {
          purityPloidy[[paste0(colnames(tumorCov)[sampleName2], "-", colnames(normalCov)[i])]] = predictWholeGenomeEvent(bAlleleFreqs[[1]], bAlleleFreqs[[2]], 
                                                                                                                         matrixOfCoverages, allowedChromsBafForThatSamples, bedFile)
        }
      } 
    }
  }
  return(purityPloidy)
}

estimatePurity <- function(delOrDup, meanValue) {
  if (delOrDup == "del") {
    purity = 4 * (meanValue) / (1 + 2 * meanValue)
  } else {
    purity = 4 * meanValue / (1 - 2 * meanValue)
  }
  return(purity)
}

