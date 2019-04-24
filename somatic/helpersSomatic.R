
findSDsOfSamples <- function(pairs, normalCov, tumorCov, bedFileForCalc, bordersOfChroms, genderOfSamplesLocal) {
  females = genderOfSamplesLocal[which(genderOfSamplesLocal == "F")]
  males = genderOfSamplesLocal[which(genderOfSamplesLocal == "M")]
  whichAreMales = which(colnames(normalCov) %in% names(males))
  whichAreFemales = which(colnames(normalCov) %in% names(females))
  chrX = which(bedFileForCalc[,1] == "chrX")
  for (i in 1:length(chrX)) {
    normalCov[chrX[i],whichAreMales] = sample(normalCov[chrX[i],whichAreFemales], length(whichAreMales), replace = T)
  }
  chrY = which(bedFileForCalc[,1] == "chrY")
  for (i in 1:length(chrY)) {
    normalCov[chrY[i],whichAreFemales] = sample(normalCov[chrY[i],whichAreMales], length(whichAreFemales), replace = T)
  }
  
  normalCovForVarianceMedians = apply(log2(normalCov), 1, median)
  
  normalCenteredAroundZero = sweep(log2(normalCov), 1, normalCovForVarianceMedians)
  normalCovForVarianceSD = apply(normalCenteredAroundZero, 2, Qn)
  normalCenteredAroundZeroAndNormalisedBySD = sweep(normalCenteredAroundZero, 2, normalCovForVarianceSD, FUN="/")
  probeVariance = apply(normalCenteredAroundZeroAndNormalisedBySD, 1, Qn)
  
  matrixOfPairs <- matrix(0, nrow=nrow(normalCov), ncol=0)
  counter = 0
  for (i in 1:ncol(normalCov)) {
    sampleNames1 <- which(pairs[,2] == colnames(normalCov)[i])
    for (sampleName1 in sampleNames1) {
      sampleName2 <- which(colnames(tumorCov) == pairs[sampleName1,1])
      new_name <- (paste(colnames(tumorCov)[sampleName2], "-",  colnames(normalCov)[i], sep=""))
      if (length(sampleName2) > 0) {
        counter = counter + 1
        matrixOfPairs <- cbind(matrixOfPairs, log2(tumorCov[,sampleName2]), log2(normalCov[,i]))
        colnames(matrixOfLogFold)[ncol(matrixOfLogFold)] <- new_name
      }
    }
  }
  sdsTum = rep(0, (ncol(matrixOfPairs) / 2))
  sdsNorm = rep(0, (ncol(matrixOfPairs) / 2))
  covNormTum = rep(0, (ncol(matrixOfPairs) / 2))
  resSDofPair = rep(0, (ncol(matrixOfPairs) / 2))
  for (j in seq(from=1, to=ncol(matrixOfPairs), by=2 )) {
    tumorS = matrixOfPairs[which(!bedFileForCalc[,1] %in% c("chrX", "chrY")),j]
    normalS = matrixOfPairs[which(!bedFileForCalc[,1] %in% c("chrX", "chrY")),j + 1]
    sdNorm = Qn(normalS)
    sdTum = Qn(tumorS)
    sdsS <- rep(0, length(bordersOfChroms) - 1)
    sdsN <- rep(0, length(bordersOfChroms) - 1)
    covNS <- rep(0, length(bordersOfChroms) - 1)
    if (length(bordersOfChroms) == 1) {
      sdsS = Qn(tumorS)
      covNS = 0
      sdsN = Qn(normalS)
    } else {
      for (i in 2:length(bordersOfChroms)) {
        if (bordersOfChroms[i] < length(tumorS)) {
          valuesBetweenBordersTum <- tumorS[bordersOfChroms[i-1]:bordersOfChroms[i]]
          valuesBetweenBordersNorm <- normalS[bordersOfChroms[i-1]:bordersOfChroms[i]]
          valuesBetweenBordersSD <- probeVariance[bordersOfChroms[i-1]:bordersOfChroms[i]]
          valuesBetweenBordersTum = (valuesBetweenBordersTum - median(valuesBetweenBordersTum)) / valuesBetweenBordersSD
          valuesBetweenBordersNorm = (valuesBetweenBordersNorm - median(valuesBetweenBordersNorm)) / valuesBetweenBordersSD
          sdsS[i - 1] = Qn(valuesBetweenBordersTum)**2
          sdsN[i - 1] = Qn(valuesBetweenBordersNorm)**2
          if (IQR(valuesBetweenBordersTum) > 0 & IQR(valuesBetweenBordersNorm) > 0) {
            covNS[i - 1] = robust_correlation_short(sqrt(sdsS[i - 1]), sqrt(sdsN[i - 1]), valuesBetweenBordersTum, valuesBetweenBordersNorm, Qn) * sqrt(sdsS[i - 1]) * sqrt(sdsN[i - 1])
          } else {
            print(i)
            covNS[i - 1] = 0
          }
        }
      }
    }
    resSDarray = (sdsS + sdsN - 2 * covNS)[(sdsS > 0)]
    indexOfMedian <- which.min(abs(resSDarray - median(resSDarray)))
    sdTum = (sdsS[sdsS > 0][indexOfMedian])
    sdNorm = (sdsN[sdsS > 0][indexOfMedian])
    covNS = (covNS[sdsS > 0][indexOfMedian])
    resSD = sqrt(resSDarray[indexOfMedian])
    sdsTum[(j + 1)/2] = sdTum
    sdsNorm[(j + 1)/2] = sdNorm
    covNormTum[(j + 1)/2] = covNS
    resSDofPair[(j + 1)/2] = resSD
  }
  return(list(rbind(sdsTum, sdsNorm, covNormTum, resSDofPair), probeVariance))
}




formilngLogFoldChange <- function(pairs, normalCov, tumorCov, currentBedFile, genderOfSamplesInCluster) {
  matrixOfLogFold <- matrix(0, nrow=nrow(normalCov), ncol=0)
  matrixOfLogFold <- foreach(i=1:ncol(normalCov), .combine="cbind") %dopar% {
    matrixOfLogFoldTmp <- matrix(0, nrow=nrow(normalCov), ncol=0)
    sampleNames1 <- which(pairs[,2] == colnames(normalCov)[i])
    for (sampleName1 in sampleNames1) {
      sampleName2 <- which(colnames(tumorCov) == pairs[sampleName1,1])
      new_name <- (paste(colnames(tumorCov)[sampleName2], "-",  colnames(normalCov)[i], sep=""))
      if (length(sampleName2) > 0) {
        matrixOfLogFoldTmp <- cbind(matrixOfLogFoldTmp, matrix(log2(tumorCov[,sampleName2]/normalCov[,i]), nrow=nrow(normalCov), ncol=1))
      }
    }
    if (ncol(matrixOfLogFoldTmp) > 0) {
      matrixOfLogFoldTmp
    } else {
      NULL
    }
  }
  counter = 1
  colnamesForMatrix <- rep(0, ncol(matrixOfLogFold))
  gendersInFormedMatrix = c()
  for (i in 1:ncol(normalCov)) {
    sampleNames1 <- which(pairs[,2] == colnames(normalCov)[i])
    for (sampleName1 in sampleNames1) {
      sampleName2 <- which(colnames(tumorCov) == pairs[sampleName1,1])
      new_name <- (paste(colnames(tumorCov)[sampleName2], "-",  colnames(normalCov)[i], sep=""))
      if (length(sampleName2) > 0) {
        colnamesForMatrix[counter] <- new_name
        gendersInFormedMatrix <- c(gendersInFormedMatrix, (genderOfSamplesInCluster[which(names(genderOfSamplesInCluster) == colnames(normalCov)[i])]))
        counter = counter + 1
      }
    }
  }
  colnames(matrixOfLogFold) = colnamesForMatrix
  
  ### THIS IS TOO SLOW! HAS TO BE RE-DONE
  uniqueChroms = unique(currentBedFile[,1])
  shifts = rep(0, nrow(matrixOfLogFold))
  for (i in 1:length(uniqueChroms)) {
    whichChrom = which(currentBedFile[,1] == uniqueChroms[i])
    matrixOfLogFoldToCheck = matrixOfLogFold[whichChrom,,drop=F]
    if (uniqueChroms[i] == "chrX") {
      if (length(which(gendersInFormedMatrix == "F")) > 10)
        matrixOfLogFoldToCheck = matrixOfLogFoldToCheck[,which(gendersInFormedMatrix == "F"),drop=F]
      else next
    }
    if (uniqueChroms[i] == "chrY") {
      if (length(which(gendersInFormedMatrix == "M")) > 10)
        matrixOfLogFoldToCheck = matrixOfLogFoldToCheck[,which(gendersInFormedMatrix == "M"),drop=F]
      else next
    }
    shiftsChrom <- apply(matrixOfLogFoldToCheck, 1, EstimateModeForNormalization)
    predictions <- runmed(shiftsChrom, k = 51)
    shifts[whichChrom] = predictions
  }
  png(paste0(opt$out, "/plot_with_shifts.png"), width=2000, height=1000)
  plot(shifts)
  lines(shifts, col="red", lwd=3)
  dev.off()
  #matrixOfLogFold <- sweep(matrixOfLogFold, 1, shifts)
  return(list(matrixOfLogFold, gendersInFormedMatrix))
}


determineSDsOfSomaticSample <- function(x, bedFile="") {
  if (!bedFile == "")
    lengtBed = bedFile[,3] - bedFile[,2]
  sdsS <- rep(0, length(bordersOfChroms) - 1)
  if (length(bordersOfChroms) == 1) {
    sdsS[1] = Qn(x)
  }
  else {
    for (i in 2:length(bordersOfChroms)) {
      valuesBetweenBorders <- x[bordersOfChroms[i-1]:bordersOfChroms[i]]
      if (!bedFile == "")
        if (unique(bedFile[bordersOfChroms[i-1]:bordersOfChroms[i],1]) %in% c("chrX","chrY")) next
      regionLengthsBetweenBorders = lengtBed[bordersOfChroms[i-1]:bordersOfChroms[i]]
      if (!bedFile == "")
        valuesBetweenBorders = valuesBetweenBorders[which(regionLengthsBetweenBorders > min(200, median(lengtBed)))]
      if (length(valuesBetweenBorders) > 1)
        sdsS[i] = Qn(valuesBetweenBorders)
    }
  }
  return(median(sdsS[which(sdsS > 0)]))
}

determineSDsOfSomaticSampleWithAllowdChroms <- function(x) {
  return(Qn(x))
}


determineSDsOfSomaticProbe <- function(x, i) {
  if (bedFile[i,1] %in% c("chrX", "chrY")) {
    x = x[which(x > median(x))]
  }
  return(Qn(x))
}


form_matrix_of_likeliks_one_sample <- function(i, vector_of_values, sds, cn_states) {
  
  vector_of_states <- cn_states
  
  sdsTmp = sds
  matrix_of_BFs = sapply(1:length(vector_of_states), function(l) {
    sds = sdsTmp
    value = return_likelik((vector_of_values - vector_of_states[l]) / (sds ) ) / (sds ) + 10^-100
    return(-2 * log(value))
  })
  return(matrix_of_BFs)
}










return_likelik <- function(x) {
  x = as.vector(x)
  x = round(abs(x * 1000)) + 1
  x = replace(x, which(x >= length(vect_of_t_likeliks)), length(vect_of_t_likeliks) - 1)
  return(vect_of_t_likeliks[x])
}

EstimateModeForNormalization <- function(x) {
  tmpx = x[which(x > -0.25 & x < 0.25)]
  if (length(tmpx) < 10) {
    return(0)
  }
  density_of_x <-  density(tmpx, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  mu
}

EstimateModeComplex <- function(x) {
  tmpx = x[which(x > -0.5 & x < 0.5)]
  if (length(tmpx) < 10) {
    tmpx = x
  }
  values <- apply(combn(tmpx, 2), 2, mean)
  density_of_x <-  density(values, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  mu
}



returnSdsForSampleAndProbe <- function(i, j) {
  sdToReturn = sdsOfSomaticSamples[j]
  return(sdToReturn * esimtatedVarianceFromSampleNoise[i] * multiplicator)
}



qcControl <- function(sam_no, toyMatrixOfLogFold, toyLocalSds, toyMultipliersDueToLog, found_CNVs, percentage) {
  copyOftoyMatrixOfLogFold <- toyMatrixOfLogFold
  pointsThatAreCNVs <- c()
  if (nrow(found_CNVs) > 0) {
    for (i in 1:nrow(found_CNVs)) {
      pointsThatAreCNVs <- c(pointsThatAreCNVs, found_CNVs[i,2]:found_CNVs[i,3])
    }
    if (length(pointsThatAreCNVs) >= nrow(copyOftoyMatrixOfLogFold) - 10) {
      return(-1)
    }
    copyOftoyMatrixOfLogFold <- copyOftoyMatrixOfLogFold[-pointsThatAreCNVs,,drop=F]
    toyLocalSds <- toyLocalSds[-pointsThatAreCNVs]
  }
  
  if (nrow(copyOftoyMatrixOfLogFold) > 0) {
    samLogFold <- copyOftoyMatrixOfLogFold[, sam_no]
    samLogFold = (samLogFold - median(samLogFold)) / toyLocalSds
    samLogFoldThreshold = quantile(samLogFold ** 2, percentage)
    samLogFold = samLogFold[which((samLogFold ** 2) < samLogFoldThreshold)]
    finalQChisq <- pchisq(sum(samLogFold ** 2), df = length(samLogFold))
    return(finalQChisq)
  } else {
    return(-1)
  }
}


returnMultiplierDueToLog <- function(cnNorm, cnTum, sdNorm, sdTum, covNT) {
  cnNorm = 2
  cnTum = 2
  resSd = (cnNorm/2)**2 * sdNorm + (cnTum/2)**2 * sdTum - 2 * (cnTum/2) * (cnNorm/2) * covNT
  if (cnTum < 0.2) {
    resSd = 2 * resSd
  }
  return(resSd)
}


returnListOfCNVsThatDoNotPass = function(found_CNVs, bafDeviationsForComparison, multiplierOfSNVsDueToMapping, bafNormalChr, bafTumorChr, 
                                         clonalityForChecking, puritiesOfStates, relativeCNumbersOfStates, bedFileForMapping, 
                                         overdispersionNormalChr, overdispersionTumorChr,
                                         pvalueShift,
                                         toyLogFoldChange,
                                         sdOfSomaticOn,
                                         sdOfSomaticOff) {
  ### CHECK BAFS HERE FOR LOW CLONALITY
  cnvsThatShowNoBAFdeviation = c()
  for (q in 1:nrow(found_CNVs)) {
    if (relativeCNumbersOfStates[found_CNVs[q,4]] < 0.5) {
      toyLogFoldChangeBetween = toyLogFoldChange[(found_CNVs[q,2] + 1):(found_CNVs[q,3] - 1)]
      if(length(which(toyLogFoldChangeBetween > log(3/4))) > 0.25 * length(toyLogFoldChangeBetween)) {
        cnvsThatShowNoBAFdeviation = c(cnvsThatShowNoBAFdeviation, q)
        print(paste("We remove potential homozygous deletion", paste0(bedFileForMapping[1,1], ":", bedFileForMapping[found_CNVs[q,2],2], "-", bedFileForMapping[found_CNVs[q,3],3]), "due to large amount of values which are too big in the middle"))
        next
      } else {
        next
      }
      
    }
    startOfCNV = as.numeric(bedFileForMapping[found_CNVs[q,2],2])
    endOfCNV <- as.numeric(bedFileForMapping[found_CNVs[q,3],3])
    coverageInsideOff = c()
    if (sdOfSomaticOff > 0) {
      coverageInsideOff = toyLogFoldChange[which(as.numeric(bedFileForMapping[,2]) >= startOfCNV & as.numeric(bedFileForMapping[,3]) <= endOfCNV & 
                                                   (as.numeric(bedFileForMapping[,3] - as.numeric(bedFileForMapping[,2]) > 10000))
      )]
    }
    coverageInsideOn = toyLogFoldChange[which(as.numeric(bedFileForMapping[,2]) >= startOfCNV & as.numeric(bedFileForMapping[,3]) <= endOfCNV & 
                                                (as.numeric(bedFileForMapping[,3] - as.numeric(bedFileForMapping[,2]) < 10000))
    )]
    trimmedCoverageInsideOn = NULL
    trimmedCoverageInsideOff = NULL
    if (length(coverageInsideOn) > 4) {
      trimmedCoverageInsideOn = trimValues(coverageInsideOn, 0.05)
    }
    if (length(coverageInsideOff) > 4) {
      trimmedCoverageInsideOff = trimValues(coverageInsideOff, 0.05)
    }
    if (length(coverageInsideOff) < 4 & length(coverageInsideOn) > 4) {
      sdOn = sd(trimmedCoverageInsideOn)
      if (sdOn > 3 * sdOfSomaticOn) {
        cnvsThatShowNoBAFdeviation = c(cnvsThatShowNoBAFdeviation, q)
        print(paste("We remove CNV", paste0(bedFileForMapping[1,1], ":", bedFileForMapping[found_CNVs[q,2],2], "-", bedFileForMapping[found_CNVs[q,3],3]),  "; level of noise", print(sdOn / sdOfSomaticOn), "due to large amount of noise in on target reads"))
        #next
      }
    } else {
      if (length(coverageInsideOff) > length(coverageInsideOn) & !is.null(trimmedCoverageInsideOff)) {
        sdOff = sd(trimmedCoverageInsideOff)
        if (sdOff > 3 * sdOfSomaticOff) {
          cnvsThatShowNoBAFdeviation = c(cnvsThatShowNoBAFdeviation, q)
          print(paste("We remove CNV", paste0(bedFileForMapping[1,1], ":", bedFileForMapping[found_CNVs[q,2],2], "-", bedFileForMapping[found_CNVs[q,3],3]), "; level of noise", print(sdOff / sdOfSomaticOff), "due to large amount of noise in off target reads"))
          #next
        }
      } else {
        if (!is.null(trimmedCoverageInsideOn)) {
          sdOn = sd(trimmedCoverageInsideOn)
          if (sdOn > 3 * sdOfSomaticOn) {
            cnvsThatShowNoBAFdeviation = c(cnvsThatShowNoBAFdeviation, q)
            print(paste("We remove CNV", paste0(bedFileForMapping[1,1], ":", bedFileForMapping[found_CNVs[q,2],2], "-", bedFileForMapping[found_CNVs[q,3],3]),"; level of noise", print(sdOn / sdOfSomaticOn), "due to large amount of noise in on target reads"))
            #next
          }
        }
      }
    }
    varsInside = which(as.numeric(bafNormalChr[,2]) >= startOfCNV & as.numeric(bafNormalChr[,3]) <= endOfCNV)
    if (as.numeric(puritiesOfStates[found_CNVs[q,4]]) <= clonalityForChecking) {
      if (length(varsInside) < 2) {
        print(paste("We remove potential CNV", paste0(bedFileForMapping[1,1], ":", bedFileForMapping[found_CNVs[q,2],2], "-", bedFileForMapping[found_CNVs[q,3],3]), "due to absence of BAF there"))
        cnvsThatShowNoBAFdeviation = c(cnvsThatShowNoBAFdeviation, q)
      } else {
        pvalsOfVariants <- rep(1, length(varsInside))
        deviation = rep(0, length(varsInside))
        for (l in 1:length(varsInside)) {
          var = varsInside[l]
          numOne = round(as.numeric(bafNormalChr[var,5]) * as.numeric(bafNormalChr[var,6]))
          numTwo = round(as.numeric(bafTumorChr[var,5]) * as.numeric(bafTumorChr[var,6]))
          refOne = as.numeric(bafNormalChr[var,6]) - numOne
          refTwo = as.numeric(bafTumorChr[var,6]) - numTwo
          overdispNorm = overdispersionNormalChr[var]
          overdispTumo = overdispersionTumorChr[var]
          deviation[l] = (abs(numTwo - multiplierOfSNVsDueToMapping * as.numeric(bafTumorChr[var,6])) / 
                            sqrt(as.numeric(bafTumorChr[var,6]) * overdispTumo * (1 - multiplierOfSNVsDueToMapping) * multiplierOfSNVsDueToMapping))
          pvalsOfVariants[l] = min(1, passPropTestVarCorrection(numOne, numTwo, refOne, refTwo, overdispNorm, overdispTumo))
        }
        wilcox.pval = wilcox.test(deviation, bafDeviationsForComparison)$p.value
        mergedPvals = pchisq((sum(log(pvalsOfVariants))*-2), df=length(pvalsOfVariants)*2, lower.tail=F)
        boxplot(deviation, bafDeviationsForComparison, main=paste(startOfCNV, endOfCNV, round(wilcox.pval, digits=4), round(mergedPvals, 4)))
        #if ((pbinom(length(which(pvalsOfVariants < 0.01)),  length(varsInside), pvalueShift, lower.tail = F) > 10 ** -4 | length(which(pvalsOfVariants < 0.01)) / length(varsInside) < 0.05) | 
        #    mergedPvals > 10 ** -4) {
        if (wilcox.pval > 0.001 | mergedPvals > 0.05) {
          cnvsThatShowNoBAFdeviation = c(cnvsThatShowNoBAFdeviation, q)
          print(paste("We remove CNV", paste0(bedFileForMapping[1,1], ":", bedFileForMapping[found_CNVs[q,2],2], "-", bedFileForMapping[found_CNVs[q,3],3]), "potential purity", puritiesOfStates[found_CNVs[q,4]], "due to 1) low clonality AND 2) absence of clear signal from BAF (p-value:", round(wilcox.pval, 4), ")"))
        }
      }
    }
    if (as.numeric(puritiesOfStates[found_CNVs[q,4]]) > clonalityForChecking) {
      if (length(varsInside) > 1) {
        pvalsOfVariants <- rep(1, length(varsInside))
        deviation = rep(0, length(varsInside))
        for (l in 1:length(varsInside)) {
          var = varsInside[l]
          numOne = round(as.numeric(bafNormalChr[var,5]) * as.numeric(bafNormalChr[var,6]))
          numTwo = round(as.numeric(bafTumorChr[var,5]) * as.numeric(bafTumorChr[var,6]))
          refOne = as.numeric(bafNormalChr[var,6]) - numOne
          refTwo = as.numeric(bafTumorChr[var,6]) - numTwo
          overdispNorm = overdispersionNormalChr[var]
          overdispTumo = overdispersionTumorChr[var]
          deviation[l] = (abs(numTwo - multiplierOfSNVsDueToMapping * as.numeric(bafTumorChr[var,6])) / 
                            sqrt(as.numeric(bafTumorChr[var,6]) * overdispTumo * (1 - multiplierOfSNVsDueToMapping) * multiplierOfSNVsDueToMapping))
          pvalsOfVariants[l] = min(1, passPropTestVarCorrection(numOne, numTwo, refOne, refTwo, overdispNorm, overdispTumo))
        }
        wilcox.pval = wilcox.test(deviation, bafDeviationsForComparison)$p.value
        mergedPvals = pchisq((sum(log(pvalsOfVariants))*-2), df=length(pvalsOfVariants)*2, lower.tail=F)
        #if (pbinom(length(which(pvalsOfVariants < 0.01)),  length(varsInside), pvalueShift, lower.tail = F) < 10 ** -4 & 
        #   length(which(pvalsOfVariants < 0.01)) / (length(varsInside)) > 0.05 &
        #   mergedPvals < 10 ** -4) {
        if (wilcox.pval < 0.001 | mergedPvals < 0.05) {
          if (q %in% cnvsThatShowNoBAFdeviation) {
            print(paste("We remain CNV", paste0(bedFileForMapping[1,1], ":", bedFileForMapping[found_CNVs[q,2],2], "-", bedFileForMapping[found_CNVs[q,3],3]), "potential purity", puritiesOfStates[found_CNVs[q,4]], " - it was filtered out but BAF shows that something is wrong (p-value:", round(wilcox.pval, 4), ")"))
            cnvsThatShowNoBAFdeviation = setdiff(cnvsThatShowNoBAFdeviation, q)
          }
        }
      }
    }
  }
  cnvsThatShowNoBAFdeviation = unique(cnvsThatShowNoBAFdeviation)
  return(cnvsThatShowNoBAFdeviation)
}


makeBarplot <- function(allPotentialPurities, found_CNVs_total, sample_name) {
  allPotentialPurities = allPotentialPurities[allPotentialPurities > 0]
  datasetForBarplot = matrix(0, nrow=4, ncol=length(unique(allPotentialPurities)))
  colnames(datasetForBarplot) = sort(unique(allPotentialPurities))
  rownames(datasetForBarplot) = c("CNeutral", "DUP", "HIGHDUP", "DEL")
  datasetForBarplotNumber = matrix(0, nrow=4, ncol=length(unique(allPotentialPurities)))
  actualCopyNumbers = as.numeric(found_CNVs_total[,4]) + as.numeric(found_CNVs_total[,5])
  colnames(datasetForBarplotNumber) = sort(unique(allPotentialPurities))
  rownames(datasetForBarplotNumber) = c("CNeutral", "DUP", "HIGHDUP", "DEL")
  for (z in 1:nrow(found_CNVs_total)) {
    dupOrDel = 1
    if (actualCopyNumbers[z] > 2 & actualCopyNumbers[z] < 5) {
      dupOrDel = 2
    } else if (actualCopyNumbers[z] > 4) {
      dupOrDel = 3
    } else if (actualCopyNumbers[z] < 2) {
      dupOrDel = 4
    }
    datasetForBarplot[dupOrDel, which(colnames(datasetForBarplot) == found_CNVs_total[z,6])] = datasetForBarplot[dupOrDel, which(colnames(datasetForBarplot) == found_CNVs_total[z,6])] + 
      as.numeric(found_CNVs_total[z,3]) - as.numeric(found_CNVs_total[z,2])
    datasetForBarplotNumber[dupOrDel, which(colnames(datasetForBarplotNumber) == found_CNVs_total[z,6])] = datasetForBarplotNumber[dupOrDel, which(colnames(datasetForBarplotNumber) == found_CNVs_total[z,6])] + 1
    
  }
  datasetForBarplot = (datasetForBarplot / 10**6)
  maxheight = max(datasetForBarplot)
  png(paste0(sample_name, "_clonalityBarplot.png"), width=2400, height=640)
  bp <- barplot(datasetForBarplot, col=c("brown","blue","darkblue","red") ,  font.axis=2, beside=T, main=paste("Presence of clones in tumor", sample_name, ", estimated purity: ", max(as.numeric(found_CNVs_total[,6]))), ylim=c(0, 1.05 * maxheight), xlab="Subclones investigated", ylab="Length, MB")
  for (z in 1:ncol(datasetForBarplotNumber)) {
    for (v in 1:nrow(datasetForBarplotNumber)) {
      if (datasetForBarplotNumber[v,z] > 0) {
        text(bp[v,z], datasetForBarplot[v,z] + 0.02 * maxheight, datasetForBarplotNumber[v,z])
      }
    }
    currentPurity = colnames(datasetForBarplot)[z]
    found_CNVs_total_LOH = found_CNVs_total[which(actualCopyNumbers== 2 & found_CNVs_total[,6] == currentPurity),,drop=F]
    if (nrow(found_CNVs_total_LOH) > 1) {
      linesToDepict = cumsum(sort(as.numeric(found_CNVs_total_LOH[, 3]) 
                                  - as.numeric(found_CNVs_total_LOH[, 2]), decreasing = T) / 10**6)[1:(nrow(found_CNVs_total_LOH) - 1)]
      for (height in linesToDepict)
        segments(bp[1,z] - 0.4, height, bp[1,z] + 0.4, height, col="white", lwd=3)
    }
    
    found_CNVs_total_dup = found_CNVs_total[which(actualCopyNumbers > 2 & as.numeric(found_CNVs_total[,4]) < 5 & found_CNVs_total[,6] == currentPurity),,drop=F]
    if (nrow(found_CNVs_total_dup) > 1) {
      linesToDepict = cumsum(sort(as.numeric(found_CNVs_total_dup[, 3]) 
                                  - as.numeric(found_CNVs_total_dup[, 2]), decreasing = T) / 10**6)[1:(nrow(found_CNVs_total_dup) - 1)]
      for (height in linesToDepict)
        segments(bp[2,z] - 0.4, height, bp[2,z] + 0.4, height, col="white", lwd=3)
    }
    
    found_CNVs_total_high_dup = found_CNVs_total[which(actualCopyNumbers > 4 & found_CNVs_total[,6] == currentPurity),,drop=F]
    if (nrow(found_CNVs_total_high_dup) > 1) {
      linesToDepict = cumsum(sort(as.numeric(found_CNVs_total_high_dup[, 3]) 
                                  - as.numeric(found_CNVs_total_high_dup[, 2]), decreasing = T) / 10**6)[1:(nrow(found_CNVs_total_high_dup) - 1)]
      for (height in linesToDepict)
        segments(bp[3,z] - 0.4, height, bp[3,z] + 0.4, height, col="white", lwd=3)
    }
    
    found_CNVs_total_del = found_CNVs_total[which(actualCopyNumbers < 2 & found_CNVs_total[,6] == currentPurity),,drop=F]
    if (nrow(found_CNVs_total_del) > 1) {
      linesToDepict = cumsum(sort(as.numeric(found_CNVs_total_del[, 3]) 
                                  - as.numeric(found_CNVs_total_del[, 2]), decreasing = T) / 10**6)[1:(nrow(found_CNVs_total_del) - 1)]
      for (height in linesToDepict)
        segments(bp[4,z] - 0.4, height, bp[4,z] + 0.4, height, col="white", lwd=3)
    }
  }
  
  dev.off()
}

makeTransparent = function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}






plotChromosomalLevelInstabs <- function(found_CNVs_total, left_borders, right_borders, ends_of_chroms, gender, sample_name) {
  found_CNVs_total[,6] = as.numeric(found_CNVs_total[,6]) / max(as.numeric(found_CNVs_total[,6]))
  majorClone = max(as.numeric(found_CNVs_total[,6]))
  linesOnBarplot = list()
  orderOfNames = c(paste0("chr", 1:22), "chrX", "chrY")
  orderInLists = c()
  for (l in 1:length(orderOfNames)) {
    orderInLists = c(orderInLists, which(names(left_borders) == orderOfNames[l]))
  }
  for (l in orderInLists) {
    startOfChr = 0
    endOfLeftArm = left_borders[[l]]
    startOfRightArm = right_borders[[l]]
    endOfRightArm = ends_of_chroms[[l]]
    nameOfChrom = names(left_borders)[l]
    linesOnBarplot[[as.character(l)]] = c(startOfChr, endOfLeftArm, startOfRightArm, endOfRightArm, nameOfChrom)
  }
  colForMajor=c("brown","blue","purple3","red")
  colForMinor = c("brown1", "lightsteelblue3", "mediumpurple1", "lightpink2")
  
  multiplicator = 80
  offsetOfSecondChr = (multiplicator / 2.5)
  widthOfLine = c(((2.3 / 20) * multiplicator), ((1.9 / 20) * multiplicator))
  pdf(file=paste0(sample_name, "_chromPlot.pdf"), width=16, height=14)
  #par(mfrow=c(2,1), mar=c(1.5, 0, 2, 1.5))
  colOfChr = c("black", "whitesmoke")
  par( mar=c(1.5, 2, 2, 1.5))
  
  chromsToAnalyse = 1:24
  
  plot(0,0, ylim=c(multiplicator - offsetOfSecondChr, multiplicator *24), xlim=c(0, max(unlist(ends_of_chroms))), col="white", xaxt="n", bty="n", axes=F, xlab="", ylab="", main=ifelse(l==1, sample_name, ""))
  
  legend("right", legend=c( "Major clone Dup 3CN", "Major clone Dup >= 4","Major clone Del",
                            "Minor clone Dup 3CN","Minor clone Dup >= 4","Minor clone Del"),
         col=c(colForMajor[2:4],colForMinor[2:4]), cex=1.8, lwd=widthOfLine, box.lty=0, 
         title=paste("Clonal fraction:", paste(round(sort(unique(as.numeric(found_CNVs_total[,6]))), digits=2),collapse="; "))
         )
  
  
  text(y = multiplicator *1:24 + offsetOfSecondChr / 2, x = rep(0 ** 7, 24), labels=orderOfNames[sort(chromsToAnalyse, decreasing = T)], pos=2, offset = 0.5)
  for (z in sort(chromsToAnalyse, decreasing = T)) {
    i = 25 - which(chromsToAnalyse == z)
    chromStructure = linesOnBarplot[[(z)]]
    for (plotNumber in 1:2) {
      if (!chromStructure[[5]] %in% c("chrX", "chrY") | (chromStructure[[5]] == "chrX" & gender == "F")) {
        segments( 0, multiplicator* i, as.numeric(chromStructure[2]), multiplicator * i,  lwd=widthOfLine[plotNumber], col=colOfChr[plotNumber])
        segments(as.numeric(chromStructure[3]),  multiplicator* i, as.numeric(chromStructure[4]), multiplicator * i,lwd=widthOfLine[plotNumber], col=colOfChr[plotNumber])
        segments(0,  multiplicator* i + offsetOfSecondChr,  as.numeric(chromStructure[2]) , multiplicator * i + offsetOfSecondChr, lwd=widthOfLine[plotNumber], col=colOfChr[plotNumber])
        segments(as.numeric(chromStructure[3]),  multiplicator* i + offsetOfSecondChr,  as.numeric(chromStructure[4]), multiplicator * i + offsetOfSecondChr,lwd=widthOfLine[plotNumber], col=colOfChr[plotNumber])
      } else {
        if ((chromStructure[[5]] == "chrX" | chromStructure[[5]] == "chrY") & gender == "M") {
          segments( 0, multiplicator* i, as.numeric(chromStructure[2]), multiplicator * i,  lwd=widthOfLine[plotNumber], col=colOfChr[plotNumber])
          segments(as.numeric(chromStructure[3]),  multiplicator* i, as.numeric(chromStructure[4]), multiplicator * i,lwd=widthOfLine[plotNumber], col=colOfChr[plotNumber])
        }
      }
    }
    # DEPICTION OF CNVs
    whichCNVsToPlot = which(found_CNVs_total[,1] == chromStructure[[5]])
    
    if (length(whichCNVsToPlot) > 0) {
      for (numOfCNV in 1:length(whichCNVsToPlot)) {
        m = whichCNVsToPlot[numOfCNV]
        particularPurity = as.numeric(found_CNVs_total[m,6])
        colorForPlotting = colForMajor
        cnvLty = 1
        cnvLwd = max(0.3, 0.9 * particularPurity) * widthOfLine[1]
        if (as.numeric(found_CNVs_total[m,6]) < majorClone - 10 ** -5) {
          colorForPlotting = colForMinor
          cnvLty = 1
        }
        cnvToPlot = found_CNVs_total[m,]
        start = as.numeric(found_CNVs_total[m,2])
        end = as.numeric(found_CNVs_total[m,3])
        copy_number_particuar_cnv_minor = as.numeric(found_CNVs_total[m,5])
        copy_number_particuar_cnv_major = as.numeric(found_CNVs_total[m,4])
        copy_number_particuar_cnv = copy_number_particuar_cnv_minor + copy_number_particuar_cnv_major
        
        colorType = c(0,0)
        if (round(copy_number_particuar_cnv_minor) < 1) {
          colorType[2] = 4
        }
        if (round(copy_number_particuar_cnv_minor) == 2) {
          colorType[2] = 2
        }
        if (round(copy_number_particuar_cnv_minor) > 2) {
          colorType[2] = 3
        }
        if (round(copy_number_particuar_cnv_major) < 1) {
          colorType[1] = 4
        }
        if (round(copy_number_particuar_cnv_major) == 2) {
          colorType[1] = 2
        }
        if (round(copy_number_particuar_cnv_major) > 2) {
          colorType[1] = 3
        }
        cnv_state = (found_CNVs_total[m,9])
        
        segments(start, multiplicator* i,  end,multiplicator * i, lwd=cnvLwd, lty = cnvLty, col=makeTransparent(colorForPlotting[colorType[1]], alpha=max(0.3, particularPurity)))
        if (colorType[2] != 0 & !(chromStructure[[5]] %in% c("chrX", "chrY") & gender == "M")) {
          segments( start, multiplicator* i + offsetOfSecondChr,  end, multiplicator * i + offsetOfSecondChr, lwd=cnvLwd, lty = cnvLty, col=makeTransparent(colorForPlotting[colorType[2]], alpha=max(0.3, particularPurity)))
        }
        
        text(y = multiplicator *i + offsetOfSecondChr / 2, x = start + (end - start) / 2, labels=paste0(copy_number_particuar_cnv), adj=c(0.5,0.5), col=ifelse(copy_number_particuar_cnv < 6, "black", "darkred"), cex=0.7)
      }
    }
    
    
    
  }
  dev.off()
}

probeLevelQC <- function(matrixOfLogFoldForCalc, sdsOfProbes, sdsOfSomaticSamples, genderOfSamplesLocal, bedFile) {
  females = which(genderOfSamplesLocal == "F")
  males = which(genderOfSamplesLocal == "M")
  #QNs <- apply(matrixOfLogFold, 1, function(x) {ifelse(length(which(x > log(1/2) & x < log(3/2))) > 10, Qn(x[which(x > log(1/2) & x < log(3/2))]), Qn(x))})
  QNs <- apply(matrixOfLogFoldForCalc, 1, Qn)
  for (i in 1:nrow(bedFile)) {
    if (bedFile[i,1] == "chrX") {
      matrixOfLogFoldForCalc[i,males] = sample(matrixOfLogFoldForCalc[i,females], length(males), replace = T)
      QNs[i] = Qn(matrixOfLogFoldForCalc[i, females])
    }
    if (bedFile[i,1] == "chrY") {
      QNs[i] = Qn(matrixOfLogFoldForCalc[i, males])
    }
  }
  ratios <- (QNs / (sdsOfProbes * median(sdsOfSomaticSamples)))
  qnRatios = Qn(ratios[which(is.finite(ratios))])
  medianRatio = median(ratios[which(is.finite(ratios))])
  plot(ratios, col=rgb(0,0,0,0.1), pch=19)
  points(ratios[which(( ratios < medianRatio - 3 * qnRatios) & !bedFile[,1] %in% c("chrY", "chrX"))] ~
           which(( ratios < medianRatio - 3 * qnRatios) & !bedFile[,1] %in% c("chrY", "chrX")), col=rgb(1,0,0,0.5), pch=19)
  sdsOfProbesCorrected = sdsOfProbes
  #sdsOfProbesCorrected[which(ratios > medianRatio + 3 * qnRatios)] = sdsOfProbesCorrected[which(ratios > medianRatio + 3 * qnRatios)] * ratios[which(ratios > medianRatio + 3 * qnRatios)]
  toRemove = which((ratios < medianRatio - 3 * qnRatios) & !bedFile[,1] %in% c("chrY", "chrX"))
  toRemove = union(toRemove, which(sdsOfProbes < 0.0001))
  return(list(toRemove, sdsOfProbesCorrected))
}



returnCoordsThatNeedToBeNull = function(bedFile, fileNameWithGermlineVars) {
  coordsToMakeNull = c()
  if(file.exists(fileNameWithGermlineVars)) {
    tryCatch({
      germlineVars <- read.table(fileNameWithGermlineVars, stringsAsFactors = F)
      if (nrow(germlineVars) > 0)
        for (j in 1:nrow(germlineVars)) {
          coordsToMakeNull = c(coordsToMakeNull, which(bedFile[,1] == germlineVars[j,1] & as.numeric(bedFile[,2]) >= as.numeric(germlineVars[j,2]) & as.numeric(bedFile[,3]) <= as.numeric(germlineVars[j,3])))
        }
    }
    , error = function(err) {
      print(err)
      
    }
    )
  }
  return(coordsToMakeNull)
}



normalizeToCommonMedian <- function(normalCov, tumorCov, bedFileForCalc, genderOfSamplesLocal) {
  normalCovForCalc = (normalCov)
  #tumorCovForCalc = (tumorCov)
  mediansOfCoverageAfterNorm = rep(0, nrow(normalCovForCalc))
  whichAutosomes = which(!bedFileForCalc[,1] %in% c("chrX", "chrY"))
  for (i in whichAutosomes) {
    mediansOfCoverageAfterNorm[i] = median(normalCovForCalc[i,])
  }
  females = genderOfSamplesLocal[which(genderOfSamplesLocal == "F")]
  males = genderOfSamplesLocal[which(genderOfSamplesLocal == "M")]
  whichAreMales = which(colnames(normalCovForCalc) %in% names(males))
  whichAreFemales = which(colnames(normalCovForCalc) %in% names(females))
  chrX = which(bedFileForCalc[,1] == "chrX")
  for (i in 1:length(chrX)) {
    normalCovForCalc[chrX[i],whichAreMales] = sample(normalCovForCalc[chrX[i],whichAreFemales], length(whichAreMales), replace = T)
  }
  for (i in chrX) {
    mediansOfCoverageAfterNorm[i] = median(normalCovForCalc[i,])
  }
  chrY = which(bedFileForCalc[,1] == "chrY")
  for (i in 1:length(chrY)) {
    normalCovForCalc[chrY[i],whichAreFemales] = sample(normalCovForCalc[chrY[i],whichAreMales], length(whichAreFemales), replace = T)
  }
  for (i in chrY) {
    mediansOfCoverageAfterNorm[i] = 2 * median(normalCovForCalc[i,])
  }
  normalCov = sweep(normalCov, 1, mediansOfCoverageAfterNorm, FUN="/")
  tumorCov = sweep(tumorCov, 1, mediansOfCoverageAfterNorm, FUN="/")
  return(list(normalCov, tumorCov))
}




findDeviationInNormalCoverage <- function(germline_sample_name, tumor_sample_name, found_CNVs_total, bedFileForCluster, tmpNormal,
                                          bedFileForClusterOff=NULL, tmpNormalOff=NULL) {
  shifts <- matrix(0, nrow=0, ncol=2)
  for (i in 1:nrow(found_CNVs_total)) {
    coordsInOn = which(bedFileForCluster[,1] == found_CNVs_total[i,1] & 
                         as.numeric(bedFileForCluster[,2]) >= as.numeric(found_CNVs_total[i,2]) & 
                         as.numeric(bedFileForCluster[,3]) <= as.numeric(found_CNVs_total[i,3]))
    valuesInSampleOn = tmpNormal[coordsInOn,which(colnames(tmpNormal) == germline_sample_name)]
    valuesCohortOn = (tmpNormal[coordsInOn,which(colnames(tmpNormal) != germline_sample_name),drop=F])
    valuesInSampleOff = c()
    valuesCohortOff = matrix(0, nrow=0, ncol=ncol(tmpNormal) - 1)
    if (!is.null(bedFileForClusterOff)) {
      coordsInOff = which(bedFileForClusterOff[,1] == found_CNVs_total[i,1] & 
                            as.numeric(bedFileForClusterOff[,2]) >= as.numeric(found_CNVs_total[i,2]) & 
                            as.numeric(bedFileForClusterOff[,3]) <= as.numeric(found_CNVs_total[i,3]))
      if (length(coordsInOff) > 0){
        valuesInSampleOff = tmpNormalOff[coordsInOff,which(colnames(tmpNormalOff) == germline_sample_name)]
        valuesCohortOff = (tmpNormalOff[coordsInOff,which(colnames(tmpNormalOff) != germline_sample_name),drop=F])
      }
    }
    valuesSample = median(c(valuesInSampleOn, valuesInSampleOff))
    valuesInSampleOff = c()
    if (!is.null(bedFileForClusterOff)) {
      if (ncol(valuesCohortOff) > 5 & nrow(valuesCohortOff) > 0) {
        valuesCohortOn = valuesCohortOn[,which(colnames(valuesCohortOn) %in% colnames(valuesCohortOff)),drop=F]
        valuesCohortOff = valuesCohortOff[,which(colnames(valuesCohortOff) %in% colnames(valuesCohortOn)),drop=F]
        valuesCohortOn = valuesCohortOn[,order(colnames(valuesCohortOn)),drop=F]
        valuesCohortOff = valuesCohortOff[,order(colnames(valuesCohortOff)),drop=F]
      } 
    } else {
      valuesCohortOff = matrix(0, nrow=0, ncol=ncol(tmpNormal) - 1)
    }
    valuesCohort <- apply(rbind(valuesCohortOn, valuesCohortOff), 2, median)
    prob = 2 * pt(-abs(   
      valuesSample - median(valuesCohort)
    ) / Qn(valuesCohort), df=length(valuesCohort)) 
    ratio = valuesSample / median(valuesCohort)
    if (found_CNVs_total[i,1] %in% c("chrX", "chrY")) ratio = 1
    if (prob < 0.01 & abs(valuesSample - median(valuesCohort)) > 0.025) {
      shifts = rbind(shifts, matrix(c(ratio, F) , ncol=2))
      print(paste0("Potential normal-specific CNVS! Shift of coverage: ", valuesSample - median(valuesCohort)))
      print(paste(found_CNVs_total[i,1:6], collapse="  "))
    } else {
      shifts = rbind(shifts, matrix(c(ratio, T) , ncol=2))
    }
  }
  return(shifts)
}


returnAreasFreeOfCNVsForAdditionalAnalysis <- function(found_CNVs_total, sample_gender, bedFileForCluster, bedFileForClusterOff=NULL) {
  uniqueChroms = unique(bedFileForCluster[,1])
  areasFreeOfCNVs <- matrix(nrow=0, ncol=5)
  colnames(areasFreeOfCNVs) = c("chr", "start", "end", "CN", "number_of_regions")
  for (i in 1:length(uniqueChroms)) {
    defaultCN = 2
    if (sample_gender == "M" & uniqueChroms[i] %in% c("chrX","chrY")) {
      defaultCN = 1
    }
    if (sample_gender == "F" & uniqueChroms[i] %in% c("chrY")) {
      next
    }
    sortedCNVs <- found_CNVs_total[which(found_CNVs_total[,1] == uniqueChroms[i]),,drop=F]
    sortedCNVs = sortedCNVs[order(as.numeric(sortedCNVs[,2])),,drop=F]
    startOfVariant = 0
    #if (nrow(sortedCNVs) > 0) {
    for (j in 1:(nrow(sortedCNVs) + 1)) {
      startOfNeutralSite = startOfVariant
      if (j <= nrow(sortedCNVs)) {
        endOfNeutralSite = as.numeric(sortedCNVs[j,2])
        endOfCNV = as.numeric(sortedCNVs[j,3])
      } else {
        endOfNeutralSite = max(as.numeric(bedFileForCluster[which(bedFileForCluster[,1] == uniqueChroms[i]),3]))
        if (!is.null(bedFileForClusterOff)) {
          endOfNeutralSite = max(endOfNeutralSite, max(as.numeric(bedFileForClusterOff[which(bedFileForClusterOff[,1] == uniqueChroms[i]),3])))
        }
      }
      formedVariant = c(uniqueChroms[i], startOfNeutralSite, endOfNeutralSite, defaultCN, length(
        which(bedFileForCluster[,1] == uniqueChroms[i] & as.numeric(bedFileForCluster[,2]) >= startOfNeutralSite & as.numeric(bedFileForCluster[,3]) <= endOfNeutralSite)))
      if (!is.null(bedFileForClusterOff)) {
        formedVariant[5] = as.numeric(formedVariant[5]) + length(
          which(bedFileForClusterOff[,1] == uniqueChroms[i] & as.numeric(bedFileForClusterOff[,2]) >= startOfNeutralSite & as.numeric(bedFileForClusterOff[,3]) <= endOfNeutralSite))
      }
      if (j <= nrow(sortedCNVs)) {
        startOfVariant = endOfCNV
      }
      areasFreeOfCNVs = rbind(areasFreeOfCNVs, formedVariant)
    }
    #}
    
  }
  areasFreeOfCNVs = areasFreeOfCNVs[which(as.numeric(areasFreeOfCNVs[,5]) > 1),,drop=F]
  return(areasFreeOfCNVs)
}

























plotFoundCNVsNew <- function(sam_no, found_CNVs, toyLogFoldChange, toyBedFile, outputFolder, chrom, cn_states, local_copy_numbers_used_major, local_copy_numbers_used_minor, purities, toySizesOfPointsFromLocalSds, plottingOfPNGs) {
  vector_of_states <- cn_states
  cnvsToOutput <- matrix(0, nrow=0, ncol=10)
  if (nrow(found_CNVs) > 0) {
    for (s in 1:nrow(found_CNVs)) {
      if(opt$debug) {
        print("Started with")
      }
      CNV_name <- paste(chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], "CN:", vector_of_states[found_CNVs[s,4]], "-2ln(loglik):", found_CNVs[s,1])
      CNV_name_to_write <- paste(colnames(toyLogFoldChange)[sam_no],  chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], "CN",vector_of_states[found_CNVs[s,4]], sep="_")
      
      vectorOfGeneNames = c()
      genesThatHasToBeSeparated = unique(toyBedFile[found_CNVs[s,2]:found_CNVs[s,3],5])
      for (i in 1:length(genesThatHasToBeSeparated)) {
        if (is.character(genesThatHasToBeSeparated[i]))
          vectorOfGeneNames = c(vectorOfGeneNames, unlist(strsplit(genesThatHasToBeSeparated[i], split=",")))
      }
      vectorOfGeneNamesTrimmed = c()
      if (length(vectorOfGeneNames) > 0) {
        for (elem in vectorOfGeneNames) {
          vectorOfGeneNamesTrimmed = c(vectorOfGeneNamesTrimmed,trimws(elem) )
        }
      }
      if (length(vectorOfGeneNamesTrimmed) > 0) {
        annotationGenes <- paste(unique(vectorOfGeneNamesTrimmed), collapse=",")
      } else {
        annotationGenes = "na"
      }
      CNVtoOut <- matrix(c(chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], 
                           local_copy_numbers_used_major[found_CNVs[s,4]], local_copy_numbers_used_minor[found_CNVs[s,4]],
                           purities[found_CNVs[s,4]],
                           vector_of_states[found_CNVs[s,4]], round(-1 * found_CNVs[s,1],0), 
                           found_CNVs[s,3] - found_CNVs[s,2] + 1, annotationGenes), nrow=1)
      if(opt$debug)
      {
        print(CNVtoOut)
      }
      cnvsToOutput = as.matrix(rbind(cnvsToOutput, CNVtoOut), ncol=10, drop=F)
      
      
      length_of_repr <- 500
      
      
      st <- found_CNVs[s,2]
      fn <- found_CNVs[s,3]
      
      pr = plottingOfPNGs
      if (pr) {
        png(filename=paste0(outputFolder, "/", paste0(CNV_name_to_write, ".png")), type="cairo",width = 1024, height = 640)
        if(opt$debug) {
          print(CNV_name_to_write)
          print(paste0(outputFolder, CNV_name_to_write))
        }
        
        
        plot(toyLogFoldChange, main=CNV_name, ylab="Copy Number", xlab=(paste("CNV within Chromosome Arm" )),
             ylim=c(-5, 5), cex=toySizesOfPointsFromLocalSds, yaxt='n')
        
        axis(side = 2, at = c(log2(cn_states/2)), labels = cn_states)
        abline(v=c(found_CNVs[s,2], found_CNVs[s,3]), col="red")
        
        
        
        abline(h=log2(c(1, 0.5, 3:10/2)),lty=2,col=c("darkgreen", rep("blue", 10)),lwd=3)
        points(found_CNVs[s,2]:found_CNVs[s,3], toyLogFoldChange[found_CNVs[s,2]:found_CNVs[s,3]],col="black", pch=21,bg=colours[found_CNVs[s,4]], cex=toySizesOfPointsFromLocalSds[found_CNVs[s,2]:found_CNVs[s,3]])
        
        ### EACH POINTS WITH DISTANCE > 10 MB WILL BE SEPARATED BY VERTICAL LINE
        distanceBetweenPoints = 10 ** 6
        for (i in 2:nrow(toyBedFile)) {
          if (toyBedFile[i,2] - toyBedFile[i - 1,2] > distanceBetweenPoints) {
            abline(v = i - 0.5, lty=2, col="grey")
          }
        }
        
        
        dev.off()
      }
      
    }
  }
  return(cnvsToOutput)
}


determine_potential_states = function(sampleLogFold, local_cn_states, sampleLogFoldOfftarget=NULL) {
  arrayOfMedians <- runmed(sampleLogFold, opt$lengthS)
  if (!is.null(sampleLogFoldOfftarget)) {
    arrayOfMedians <- c(arrayOfMedians, runmed(sampleLogFoldOfftarget, opt$lengthS))
  }
  diffsFromCoverage <- sapply(1:length(local_cn_states), function(i) {min(abs(log2(local_cn_states[i] / 2) - (arrayOfMedians)))})
  blocked_states = setdiff(which(diffsFromCoverage > log2(1.05)), c(1,2))
  return(blocked_states)
}



find_baseline_level <- function(allowedChromsBafSample, matrixOfLogFoldSample, bedFileForCluster, matrixOfLogFoldOffSample=NULL, bedFileForClusterOff=NULL) {
  if (is.null(matrixOfLogFoldOffSample)) {
    allowedChromosomesAutosomesOnly = c()
    for (allowedArm in allowedChromsBafSample) {
      splittedValue <- strsplit(allowedArm, "-")
      chrom = splittedValue[[1]][1]
      if (!chrom %in% c("chrY", "Y", "chrX", "X")) {
        startOfArm = as.numeric(splittedValue[[1]][2])
        endOfArm = as.numeric(splittedValue[[1]][3])
        allowedChromosomesAutosomesOnly = union(allowedChromosomesAutosomesOnly, which(bedFileForCluster[,1] == chrom &
                                                                                         bedFileForCluster[,2] >= startOfArm &
                                                                                         bedFileForCluster[,3] <= endOfArm))
      }
    }
    lengthOfRolling = 51
    matrixOfLogFoldAllowedChrom = matrixOfLogFoldSample[allowedChromosomesAutosomesOnly ]
    
    smoothedLogFold = runmed(matrixOfLogFoldAllowedChrom, k = lengthOfRolling)
    clusteredResult <- densityMclust(smoothedLogFold[which(smoothedLogFold > log2(2/8))])
    print("Mclust finished")
    bigClusters <- which(clusteredResult$parameters$pro > 0.25)
    if (length(bigClusters) == 0) {
      shiftOfCoverage <- median(matrixOfLogFold[allowedChromosomesAutosomesOnly])
    } else {
      shiftOfCoverage = min(clusteredResult$parameters$mean[bigClusters])
    }
  } else {
    globalBed <- rbind(bedFileForCluster, bedFileForClusterOff)
    globalLogFold <- c( matrixOfLogFoldSample, matrixOfLogFoldOffSample)
    globalLogFold = globalLogFold[order(globalBed[,1], as.numeric(globalBed[,2]))]
    globalBed = globalBed[order(globalBed[,1], as.numeric(globalBed[,2])),]
    allowedChromosomesAutosomesOnly = c()
    smoothedLogFold= c()
    for (allowedArm in allowedChromsBafSample) {
      splittedValue <- strsplit(allowedArm, "-")
      chrom = splittedValue[[1]][1]
      if (!chrom %in% c("chrY", "Y", "chrX", "X")) {
        startOfArm = as.numeric(splittedValue[[1]][2])
        endOfArm = as.numeric(splittedValue[[1]][3])
        
        allowedChromosomesAutosomesOnly = union(allowedChromosomesAutosomesOnly, which(globalBed[,1] == chrom &
                                                                                         as.numeric(globalBed[,2]) >= startOfArm &
                                                                                         as.numeric(globalBed[,3]) <= endOfArm))
        lengthOfRolling = min(21, round((length(allowedChromosomesAutosomesOnly))/5))
        runmedian <- function(k, vec) {sapply(1:length(vec), function(i) {median(vec[max(1,i-k):min(length(vec),i+k)])})}
        #smoothedLogFold = c(smoothedLogFold, runmed(globalLogFold[which(globalBed[,1] == chrom &
        #                                                                  as.numeric(globalBed[,2]) >= startOfArm &
        #                                                                  as.numeric(globalBed[,3]) <= endOfArm)], k = lengthOfRolling, endrule="constant"))
        smoothedLogFold = c(smoothedLogFold, runmedian(lengthOfRolling, globalLogFold[which(globalBed[,1] == chrom &
                                                                                              as.numeric(globalBed[,2]) >= startOfArm &
                                                                                              as.numeric(globalBed[,3]) <= endOfArm)]))
      }
    }
    #globalLogFoldAllowedChroms = globalLogFold[allowedChromosomesAutosomesOnly]
    #smoothedLogFold = runmed(globalLogFoldAllowedChroms, k = lengthOfRolling)
    clusteredResult <- densityMclust(smoothedLogFold[which(smoothedLogFold > log2(3/8))], model="E")
    bigClusters <- which(clusteredResult$parameters$pro > 0.25)
    if (length(bigClusters) == 0) {
      shiftOfCoverage <- median(globalLogFold[allowedChromosomesAutosomesOnly])
    } else {
      shiftOfCoverage = clusteredResult$parameters$mean
      if (length(bigClusters) > 0)
        for (bC in bigClusters) {
          currentLocation = shiftOfCoverage[bC]
          diffs = abs(clusteredResult$parameters$mean - currentLocation)
          shiftOfCoverage[bC] = clusteredResult$parameters$mean[which(diffs < 0.035)] * clusteredResult$parameters$pro[which(diffs < 0.035)]
        }
      shiftOfCoverage = shiftOfCoverage[bigClusters]
    }
    print(paste0("Mass of clusters for finding diploid state: ", paste(clusteredResult$parameters$pro[bigClusters], sep=";")))
  }
  return(shiftOfCoverage)
}
