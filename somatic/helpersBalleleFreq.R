
likelihoodOfSNV <- function(a, b, p, overdispersionValue) {
  
  if (is.na(overdispersionValue)) {
    value = dbinom(a, size=b, prob=p)
  } else{
    varBinom = p * (1-p) * b * overdispersionValue
    dist = (a - p * b) / sqrt(varBinom)
    value = return_likelik(dist)
  }
    
  if (is.nan(value) | value < 10**-15) return(10**-15)
  return((value))
}


passPropTest <- function(numOne, numTwo, refOne, refTwo) {
  matrixToTest <- matrix(c(numOne, refOne, numTwo, refTwo), nrow = 2)
  pval=0.0
  if (sum(matrixToTest <= 10) > 0) {
    pval=fisher.test(matrixToTest)$p.value
    return(pval)
  } else {
    
    pval=prop.test(matrixToTest)$p.value
    return(pval)
  }
}

passPropTestVarCorrection <- function(numOne, numTwo, refOne, refTwo, overdispNorm, overdispTumo) {
  pa = numOne / (refOne + numOne)
  pb = numTwo / (refTwo + numTwo)
  overallProp = (numOne + numTwo) / (refOne + refTwo + numOne + numTwo)
  commonVariance = ( (overallProp * (1 - overallProp) * overdispNorm + overallProp * (1 - overallProp) * overdispTumo) / 2 ) 
  z = (pa - pb) / sqrt(  commonVariance * ( 1 / (refOne + numOne) + 1 /  (refTwo + numTwo))  )
  
  pval = 2 * pnorm( -abs(z))
  return(pval)
}






determineHeterozygousPositions <- function(freq, depth, probAB=0.48, threshold = 0.01) {
  if (freq < probAB - 0.1 | freq > probAB + 0.1) return(F)
  prob = pbinom(round(freq * depth), prob=probAB, size=depth)
  if (is.na(prob)) {
    print("NA in B-allele")
    print(freq)
    print(depth)
    return(F)
  }
  if ((prob > threshold & prob < 1 - threshold)) {
    return(T)
  } else {
    return(F)
  }
}



determineHeterozygousPositionsOverdispersed <- function(freq, depth, probAB=0.48, overdispersionFactors) {
  if (freq < probAB - 0.1 | freq > probAB + 0.1) return(F)
  standardDeviation = sqrt(probAB * (1 - probAB) * depth * overdispersionFactors)
  prob = 2 * pnorm( -abs(round(freq * depth) - round(probAB * depth)) / standardDeviation)
  if (is.na(prob)) {
    print("NA in B-allele")
    print(freq)
    print(depth)
    return(F)
  }
  if ((prob > 0.05)) {
    return(T)
  } else {
    return(F)
  }
}

roundUpToDegree = function(num, digits) {
  num = round(num * degree) / degree
}

whichPTouse = function(purities, cnstates, statesUsed, multiplierOfSNVsDueToMapping, degreeOfRoughness) {
  multiplierDueToMapping = multiplierOfSNVsDueToMapping / 0.5
  whichPUsed = list()
  roundUpToDegree = function(num, digits) {
    round(num * digits) / digits
  }
  probs = c()
  upperThreshold = 0.99
  lowerThreshold = 0.01
  for (j in 1:length(cnstates)) {
    cn = cnstates[j]
    stateUsed = statesUsed[j]
    pur = purities[j]
    
    pUsed1 = NULL
    pUsed2 = NULL
    pUsed3 = NULL
    overallNumberOfReads <- (1 - pur) * 2 + pur * (cn)
    if (cn != 0 & cn != 2 & stateUsed == "CNV") {
      pUsed1=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, ((1 - pur) + pur * (cn - 1)) / overallNumberOfReads), digits=degreeOfRoughness)
      pUsed2=roundUpToDegree(multiplierDueToMapping * max(0.01, 1 / overallNumberOfReads), digits=degreeOfRoughness)
    } else if (cn == 0) {
      if (pur < 0.8) {
        pUsed1 = roundUpToDegree(0.5 * multiplierDueToMapping, degreeOfRoughness)
      } else {
        pUsed1 = roundUpToDegree(multiplierDueToMapping * lowerThreshold, degreeOfRoughness)
        pUsed2 = roundUpToDegree(multiplierDueToMapping * upperThreshold, degreeOfRoughness)
        pUsed3 = roundUpToDegree(0.5 * multiplierDueToMapping, degreeOfRoughness)
      }
    } else if (stateUsed %in% c("LOH", "LOHDup")) {
      if (cn == 2) {
        pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, 0.5 - pur/2), digits=degreeOfRoughness)
        pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, 0.5 + pur/2), digits=degreeOfRoughness)
      } else if (cn == 3) {
        pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, (1 + 2 * pur) / (2 + pur)), digits=degreeOfRoughness)
        pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, (1 - pur) / (2 + pur)), digits=degreeOfRoughness)
      } else if (cn == 4) {
        pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, (1 + 3 * pur) / (2 + 2 * pur)), digits=degreeOfRoughness)
        pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, (1 - pur) / (2 + 2 * pur)), digits=degreeOfRoughness)
      }
    } else if (stateUsed == "CNVboth" | stateUsed == "normal") {
      pUsed1 = roundUpToDegree((multiplierDueToMapping * 0.5), digits=degreeOfRoughness)
    } else if (stateUsed == "CNVcomplex2" | stateUsed == "CNVcomplex3") {
      if (cn == 5) {
        pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, (1 + pur) / (2 + 3 * pur)), digits=degreeOfRoughness)
        pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, (1 + 2 * pur) / (2 + 3* pur)), digits=degreeOfRoughness)
        
      } else if (cn == 6) {
        pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, (1 + 1 * pur) / (2 + 4 * pur)), digits=degreeOfRoughness)
        pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, (1 + 3 * pur) / (2 + 4 * pur)), digits=degreeOfRoughness)
      } else if (cn == 7) {
        if (stateUsed == "CNVcomplex3") {
          pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, (1 + 3 * pur) / (2 + 5 * pur)), digits=degreeOfRoughness)
          pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, (1 + 2 * pur) / (2 + 5 * pur)), digits=degreeOfRoughness)
        } else {
          pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, (1 + 4 * pur) / (2 + 5 * pur)), digits=degreeOfRoughness)
          pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, (1 + 1 * pur) / (2 + 5 * pur)), digits=degreeOfRoughness)
        }
      } else if (cn == 8) {
        if (stateUsed == "CNVcomplex3") {
          pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, (1 + 4 * pur) / (2 + 6 * pur)), digits=degreeOfRoughness)
          pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, (1 + 2 * pur) / (2 + 6 * pur)), digits=degreeOfRoughness)
        } else {
          pUsed1=roundUpToDegree(multiplierDueToMapping * max(lowerThreshold, (1 + 1 * pur) / (2 + 6 * pur)), digits=degreeOfRoughness)
          pUsed2=roundUpToDegree(multiplierDueToMapping * min(upperThreshold, (1 + 5 * pur) / (2 + 6 * pur)), digits=degreeOfRoughness)
        }
      }
    } 
    whichProbUsed = c()
    if (!is.null(pUsed1)) {
      probs <- c(probs, pUsed1)
      whichProbUsed <- c(whichProbUsed, pUsed1)
    }
    if (!is.null(pUsed2)) {
      probs <- c(probs, pUsed2)
      whichProbUsed <- c(whichProbUsed, pUsed2)
    }
    if (!is.null(pUsed3)) {
      probs <- c(probs, pUsed3)
      whichProbUsed <- c(whichProbUsed, pUsed3)
    }
    whichPUsed[[j]] = whichProbUsed
  } 
  return(whichPUsed)
}

fillInPList = function(value, depth, whichPUsed, overdispersionValue) {
  pList = list()
  probs = sort(unique(unlist(whichPUsed)))
  distances = abs(probs - value / depth)
  startingCoords = which.min(distances)
  pList[[as.character(probs[startingCoords])]] = likelihoodOfSNV(value, depth, probs[startingCoords], overdispersionValue)
  mimimumLikelik = 10**-15
  thresholdForComparison = 10**-14
  if (startingCoords < length(probs)) {
    meaningful = T
    for (i in (startingCoords + 1):length(probs)) {
      if (meaningful) {
        likelik = likelihoodOfSNV(value, depth, probs[i], overdispersionValue)
      } else {
        likelik = mimimumLikelik
      }
      pList[[as.character(probs[i])]] = likelik
      if (likelik < thresholdForComparison) {
        meaningful = F
      }
    }
  }
  if (startingCoords > 1) {
    meaningful = T
    for (i in startingCoords:1) {
      if (meaningful) {
        likelik = likelihoodOfSNV(value, depth, probs[i], overdispersionValue)
      } else {
        likelik = mimimumLikelik
      }
      pList[[as.character(probs[i])]] = likelik
      if (likelik < thresholdForComparison) {
        meaningful = F
      }
    }
  }
  return(pList)
}

likelihoodOfSNVBasedOnCN <- function(whichPUsed, pList) {
  sumOfLikeliks = 0
  for (pUsed in whichPUsed) {
    sumOfLikeliks = sumOfLikeliks + pList[[as.character(pUsed)]]
  }
  finalLikelihood = log(sumOfLikeliks / length(whichPUsed))
  return(finalLikelihood)
}




extractVariancesFromBAF <- function(bafTable, expectedValue) {
  freqs = as.numeric(bafTable[,5])
  depths = as.numeric(bafTable[,6])
  pbinoms <- sapply(1:length(freqs), function(i) {determineHeterozygousPositions(freqs[i], depths[i], expectedValue, 0.005)})
  depths = depths[which(pbinoms == T)]
  freqs = freqs[which(pbinoms == T)]
  if (length(depths) <= 100) {
    return(cbind(rep(1, length(depths)), depths))
  }
  
  clustering = densityMclust(freqs, modelNames=c("E"))
  dists <- (abs(clustering$parameters$mean - expectedValue))
  sortedDists <- sort(dists)
  if (length(dists) == 1) {
    clustersToUse = 1
  } else {
    for (i in 2:length(dists)) {
      clustersToUse = which(dists < sortedDists[i])
      if (length(which(clustering$classification %in% clustersToUse)) > 100 & sortedDists[i] > 0.025) {
        break
      }
    }
  }
  forEstimation = which(clustering$classification %in% clustersToUse)
  forEstimationFreqs = freqs[forEstimation]
  #plot(density(forEstimationFreqs, bw="SJ"))
  forEstimationDepths = depths[forEstimation]
  orderOfDepths = order(forEstimationDepths)
  forEstimationFreqs = forEstimationFreqs[orderOfDepths]
  forEstimationDepths = forEstimationDepths[orderOfDepths]
  if (length(forEstimationFreqs) < 100) {
    return(cbind(rep(1, length(depths)), depths))
  } else {
    vecOfSDs <- rep(0, length(forEstimationFreqs))
    rollingLength = min(30, round(length(forEstimationFreqs) / 10))
    overdispersion = rep(0, length(forEstimationFreqs))
    for (i in 1:length(vecOfSDs)) {
      if (length(which(abs(forEstimationDepths - forEstimationDepths[i]) < 0.001)) > 50) {
        vecOfSDs[i] = Qn(forEstimationFreqs[which(abs(forEstimationDepths - forEstimationDepths[i]) < 0.001)])
      } else {
        vecOfSDs[i] = Qn(forEstimationFreqs[max(1,i-rollingLength):min(length(vecOfSDs), i+rollingLength)] * forEstimationDepths[i]) ** 2
      }
      varPredictedByBinom = forEstimationDepths[i] * expectedValue * (1 - expectedValue)
      overdispersion[i] = (vecOfSDs[i]) / varPredictedByBinom
    }
    overdispersion[which(overdispersion < 1)] = 1
    overdispersion[which(overdispersion >4)] = 4
    #plot(overdispersion ~ forEstimationDepths)
    return(cbind(overdispersion, forEstimationDepths))
  }
}