
findSDsOfSamples <- function(pairs, normalCov, tumorCov, bedFileForCalc, bordersOfChroms, genderOfSamplesLocal) {
  normalCovForVarianceMedians = apply(log2(normalCov), 1, median)
  normalCenteredAroundZero = sweep(log2(normalCov), 2, normalCovForVarianceMedians)
  normalCovForVarianceSD = apply(normalCenteredAroundZero, 2, Qn)
  normalCenteredAroundZeroAndNormalisedBySD = sweep(normalCenteredAroundZero, 1, normalCovForVarianceSD, FUN="/")
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


# # CODE FOR CHECKING VALIDITY OF
# vars <- matrix(0, ncol=3, nrow=0)
# finalCohort <- matrix(0, ncol=0, nrow=nrow(matrixOfLogFold))
# for (sampleNo in 1:ncol(matrixOfLogFold)) {
#   testSample <- matrix(0, nrow=length(probeVariance), ncol=2)
#   testSample1 <- rep(0, length(probeVariance))
#   for (i in 1:length(probeVariance)) {
#     varTumNow = (sdsTum[sampleNo]) * probeVariance[i] ** 2
#     varNormNow = (sdsNorm[sampleNo]) * probeVariance[i] ** 2
#     covNow = covNormTum[sampleNo] * probeVariance[i] ** 2
#     testSample[i,] = MASS::mvrnorm(1, mu=c(0,0), Sigma = matrix(c(varTumNow, covNow, covNow, varNormNow), nrow=2))
#     testSample1[i] = rnorm(1, sd=sqrt(varTumNow + varNormNow - 2 * covNow))
#   }
#   finalCohort = cbind(finalCohort, testSample1)
#   vars <- rbind(vars, c(Qn(testSample[,1] - testSample[,2]), Qn(testSample1), Qn(matrixOfLogFold[,sampleNo])))
# }
# QnsFake <- apply(finalCohort, 1, Qn)
# QnsR <- apply(matrixOfLogFold, 1, Qn)
# summary(QnsFake / QnsR)

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
  for (i in 1:ncol(normalCov)) {
    sampleNames1 <- which(pairs[,2] == colnames(normalCov)[i])
    for (sampleName1 in sampleNames1) {
      sampleName2 <- which(colnames(tumorCov) == pairs[sampleName1,1])
      new_name <- (paste(colnames(tumorCov)[sampleName2], "-",  colnames(normalCov)[i], sep=""))
      if (length(sampleName2) > 0) {
        colnamesForMatrix[counter] <- new_name
        counter = counter + 1
      }
    }
  }
  colnames(matrixOfLogFold) = colnamesForMatrix
  
  ### THIS IS TOO SLOW! HAS TO BE RE-DONE
  uniqueChroms = unique(currentBedFile[,1])
  shifts = rep(0, nrow(matrixOfLogFold))
  for (i in 1:length(uniqueChroms)) {
    whichChrom = which(currentBedFile[,] == uniqueChroms[i])
    matrixOfLogFoldToCheck = matrixOfLogFold[whichChrom,]
    if (uniqueChroms == "chrX") {
      if (length(which(genderOfSamplesInCluster == "F")) > 10)
        matrixOfLogFoldToCheck = matrixOfLogFoldToCheck[which(genderOfSamplesInCluster == "F")]
      else next
    }
    if (uniqueChroms == "chrY") {
      if (length(which(genderOfSamplesInCluster == "M")) > 10)
        matrixOfLogFoldToCheck = matrixOfLogFoldToCheck[which(genderOfSamplesInCluster == "M")]
      else next
    }
    shiftsChrom <- apply(matrixOfLogFoldToCheck, 1, EstimateModeForNormalization)
    fit <- loess(shiftsChrom ~ c(1:length(shiftsChrom)), span=0.8)
    predictions <- predict(fit, 1:length(shiftsChrom))
    shifts[whichChrom] = predictions
  }
  shiftsAll <- apply(matrixOfLogFold, 1, EstimateModeForNormalization)
  png(paste0("plot_with_shifts.png"), width=2000, height=1000)
  plot(shiftsAll)
  lines(shifts, col="red", lwd=3)
  dev.off()
  #matrixOfLogFold <- sweep(matrixOfLogFold, 1, shifts)
  return(list(matrixOfLogFold))
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


form_matrix_of_likeliks_one_sample <- function(i, j, vector_of_values, sds, cn_states, multipliersDueToLog) {

  vector_of_states <- cn_states
  matrix_of_BFs <- matrix(0, nrow=(j - i + 1), ncol=length(vector_of_states))
  start <- 1
  end <- j - i + 1
  
  homozygousDelSD = 0.5 / 10
  sdsTmp = sds
  matrix_of_BFs = sapply(1:ncol(matrix_of_BFs), function(l) {
    sds = sdsTmp
    if (vector_of_states[l] < 0.5) {
      sds[which(sds < homozygousDelSD)] = homozygousDelSD
    }
    value = return_likelik((vector_of_values - vector_of_states[l]) / (sds * multipliersDueToLog[l]) ) / (sds * multipliersDueToLog[l]) + 10^-100
    return(-2 * log(value))
  })
  return(matrix_of_BFs)
}




plotFoundCNVs <- function(found_CNVs, toyLogFoldChange, toyBedFile, outputFolder, chrom, cn_states, copy_numbers_used, purities, local_cnv_states, toySizesOfPointsFromLocalSds, plottingOfPNGs) {
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
                           copy_numbers_used[found_CNVs[s,4]], purities[found_CNVs[s,4]],
                           vector_of_states[found_CNVs[s,4]], round(-1 * found_CNVs[s,1],0), 
                           found_CNVs[s,3] - found_CNVs[s,2] + 1,
                           local_cnv_states[found_CNVs[s,4]], annotationGenes), nrow=1)
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


returnListOfCNVsThatDoNotPass = function(foundCNVs, bafNormalChr, bafTumorChr, clonalityForChecking, puritiesOfStates, bedFileForMapping, overdispersionNormalChr, overdispersionTumorChr) {
  ### CHECK BAFS HERE FOR LOW CLONALITY
  cnvsThatShowNoBAFdeviation = c()
  for (q in 1:nrow(found_CNVs)) {
    if (puritiesOfStates[found_CNVs[q,4]] > clonalityForChecking) next
    
    startOfCNV = as.numeric(bedFileForMapping[found_CNVs[q,2],2])
    endOfCNV <- as.numeric(bedFileForMapping[found_CNVs[q,3],3])
    varsInside = which(as.numeric(bafNormalChr[,2]) >= startOfCNV & as.numeric(bafNormalChr[,3]) <= endOfCNV)
    if (length(varsInside) < 10) {
      cnvsThatShowNoBAFdeviation = c(cnvsThatShowNoBAFdeviation, q)
    } else {
      pvalsOfVariants <- rep(1, length(varsInside))
      for (l in 1:length(varsInside)) {
        var = varsInside[l]
        numOne = round(as.numeric(bafNormalChr[var,5]) * as.numeric(bafNormalChr[var,6]))
        numTwo = round(as.numeric(bafTumorChr[var,5]) * as.numeric(bafTumorChr[var,6]))
        refOne = as.numeric(bafNormalChr[var,6]) - numOne
        refTwo = as.numeric(bafTumorChr[var,6]) - numTwo
        overdispNorm = overdispersionNormalChr[var]
        overdispTumo = overdispersionTumorChr[var]
        pvalsOfVariants[l] = passPropTestVarCorrection(numOne, numTwo, refOne, refTwo, overdispNorm, overdispTumo)
      }
      if (length(which(pvalsOfVariants < 0.05)) < 0.05 * length(varsInside)) {
        cnvsThatShowNoBAFdeviation = c(cnvsThatShowNoBAFdeviation, q)
      }
    }
  }
  return(cnvsThatShowNoBAFdeviation)
}


makeBarplot <- function(allPotentialPurities, found_CNVs_total) {
  allPotentialPurities = allPotentialPurities[allPotentialPurities > 0]
  datasetForBarplot = matrix(0, nrow=4, ncol=length(unique(allPotentialPurities)))
  colnames(datasetForBarplot) = sort(unique(allPotentialPurities))
  rownames(datasetForBarplot) = c("CNeutral", "DUP", "HIGHDUP", "DEL")
  datasetForBarplotLikelik = matrix(0, nrow=4, ncol=length(unique(allPotentialPurities)))
  colnames(datasetForBarplotLikelik) = sort(unique(allPotentialPurities))
  rownames(datasetForBarplotLikelik) = c("CNeutral", "DUP", "HIGHDUP", "DEL")
  datasetForBarplotNumber = matrix(0, nrow=4, ncol=length(unique(allPotentialPurities)))
  for (z in 1:nrow(found_CNVs_total)) {
    dupOrDel = 1
    if (as.numeric(found_CNVs_total[z,4]) > 2 & as.numeric(found_CNVs_total[z,4]) < 5) {
      dupOrDel = 2
    } else if (as.numeric(found_CNVs_total[z,4]) > 4) {
      dupOrDel = 3
    } else if (as.numeric(found_CNVs_total[z,4]) < 2) {
      dupOrDel = 4
    }
    datasetForBarplot[dupOrDel, which(colnames(datasetForBarplot) == found_CNVs_total[z,5])] = datasetForBarplot[dupOrDel, which(colnames(datasetForBarplot) == found_CNVs_total[z,5])] + 
      as.numeric(found_CNVs_total[z,3]) - as.numeric(found_CNVs_total[z,2])
    datasetForBarplotLikelik[dupOrDel, which(colnames(datasetForBarplotLikelik) == found_CNVs_total[z,5])] = datasetForBarplotLikelik[dupOrDel, which(colnames(datasetForBarplotLikelik) == found_CNVs_total[z,5])] + 
      as.numeric(found_CNVs_total[z,7])
    datasetForBarplotNumber[dupOrDel, which(colnames(datasetForBarplotLikelik) == found_CNVs_total[z,5])] = datasetForBarplotNumber[dupOrDel, which(colnames(datasetForBarplotLikelik) == found_CNVs_total[z,5])] + 1
  }
  datasetForBarplot = (datasetForBarplot / 10**6)
  maxheight = max(datasetForBarplot)
  png(paste0(sample_name, "_clonalityBarplot.png"), width=2400, height=640)
  bp <- barplot(datasetForBarplot, col=c("brown","blue","darkblue","red") ,  font.axis=2, beside=T, main=paste("Presence of clones in tumor", sample_name, ", estimated purity: ", max(as.numeric(found_CNVs_total[,5]))), ylim=c(0, 1.05 * maxheight), xlab="Subclones investigated", ylab="Length, MB")
  for (z in 1:ncol(datasetForBarplotNumber)) {
    for (v in 1:nrow(datasetForBarplotNumber)) {
      if (datasetForBarplotNumber[v,z] > 0) {
        text(bp[v,z], datasetForBarplot[v,z] + 0.02 * maxheight, datasetForBarplotNumber[v,z])
      }
    }
    currentPurity = colnames(datasetForBarplot)[z]
    found_CNVs_total_LOH = found_CNVs_total[which(as.numeric(found_CNVs_total[,4]) == 2 & found_CNVs_total[,5] == currentPurity),,drop=F]
    if (nrow(found_CNVs_total_LOH) > 1) {
      linesToDepict = cumsum(sort(as.numeric(found_CNVs_total_LOH[, 3]) 
                                  - as.numeric(found_CNVs_total_LOH[, 2]), decreasing = T) / 10**6)[1:(nrow(found_CNVs_total_LOH) - 1)]
      for (height in linesToDepict)
        segments(bp[1,z] - 0.4, height, bp[1,z] + 0.4, height, col="white", lwd=3)
    }
    
    found_CNVs_total_dup = found_CNVs_total[which(as.numeric(found_CNVs_total[,4]) > 2 & as.numeric(found_CNVs_total[,4]) < 5 & found_CNVs_total[,5] == currentPurity),,drop=F]
    if (nrow(found_CNVs_total_dup) > 1) {
      linesToDepict = cumsum(sort(as.numeric(found_CNVs_total_dup[, 3]) 
                                  - as.numeric(found_CNVs_total_dup[, 2]), decreasing = T) / 10**6)[1:(nrow(found_CNVs_total_dup) - 1)]
      for (height in linesToDepict)
        segments(bp[2,z] - 0.4, height, bp[2,z] + 0.4, height, col="white", lwd=3)
    }
    
    found_CNVs_total_high_dup = found_CNVs_total[which(as.numeric(found_CNVs_total[,4]) > 4 & found_CNVs_total[,5] == currentPurity),,drop=F]
    if (nrow(found_CNVs_total_high_dup) > 1) {
      linesToDepict = cumsum(sort(as.numeric(found_CNVs_total_high_dup[, 3]) 
                                  - as.numeric(found_CNVs_total_high_dup[, 2]), decreasing = T) / 10**6)[1:(nrow(found_CNVs_total_high_dup) - 1)]
      for (height in linesToDepict)
        segments(bp[3,z] - 0.4, height, bp[3,z] + 0.4, height, col="white", lwd=3)
    }
    
    found_CNVs_total_del = found_CNVs_total[which(as.numeric(found_CNVs_total[,4]) < 2 & found_CNVs_total[,5] == currentPurity),,drop=F]
    if (nrow(found_CNVs_total_del) > 1) {
      linesToDepict = cumsum(sort(as.numeric(found_CNVs_total_del[, 3]) 
                                  - as.numeric(found_CNVs_total_del[, 2]), decreasing = T) / 10**6)[1:(nrow(found_CNVs_total_del) - 1)]
      for (height in linesToDepict)
        segments(bp[4,z] - 0.4, height, bp[4,z] + 0.4, height, col="white", lwd=3)
    }
  }
  
  dev.off()
}