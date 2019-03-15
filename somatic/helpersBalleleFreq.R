
likelihoodOfSNV <- function(a, b, p, overdispersionValue) {
  
  if (is.na(overdispersionValue)) {
    value = dbinom(a, size=b, prob=p)
  } else{
    varBinom = p * (1-p) * b * overdispersionValue
    dist = (a - p * b) / sqrt(varBinom)
    value = dnorm(dist)
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
  z = (pa - pb) / sqrt( ((overdispNorm) * overallProp) / (refOne + numOne) + ((overdispTumo) * overallProp) / (refTwo + numTwo) )
  pval = 2 * pnorm( -abs(z))
  return(pval)
}






determineHeterozygousPositions <- function(freq, depth, probAB=0.48) {
  prob = pbinom(round(freq * depth), prob=probAB, size=depth)
  if (is.na(prob)) {
    print("NA in B-allele")
    print(freq)
    print(depth)
    return(F)
  }
  if ((prob > 0.01 & prob < 0.99)) {
    return(T)
  } else {
    return(F)
  }
}

determineHeterozygousPositionsOverdispersed <- function(freq, depth, probAB=0.48, overdispersionFactors) {
  standardDeviation = sqrt(probAB * (1 - probAB) * depth * overdispersionFactors)
  prob = 2 * pnorm( -abs(round(freq * depth) - round(probAB * depth)) / standardDeviation)
  if (is.na(prob)) {
    print("NA in B-allele")
    print(freq)
    print(depth)
    return(F)
  }
  if ((prob > 0.1)) {
    return(T)
  } else {
    return(F)
  }
}




likelihoodOfSNVBasedOnCN <- function(value, depth, pur, cn, stateUsed, multiplierOfSNVsDueToMapping, pList, overdispersionValue) {
  multiplierDueToMapping = multiplierOfSNVsDueToMapping / 0.5
  overallNumberOfReads <- (1 - pur) * 2 + pur * (cn)
  pListChanged = F
  if (cn != 0 & cn != 2 & stateUsed == "CNV") {
    numberOfReadsSupportiveOne1 <- (1 - pur) + pur * (cn - 1)
    numberOfReadsSupportiveOne2 <- 1 # (1 - pur) + pur * 1
    pUsed1=round(multiplierDueToMapping * min(0.99, numberOfReadsSupportiveOne1 / overallNumberOfReads), digits=2)
    pUsed2=round(multiplierDueToMapping * max(0.01, numberOfReadsSupportiveOne2 / overallNumberOfReads), digits=2)
    if (!as.character(pUsed1) %in% names(pList)) {
      pList[[as.character(pUsed1)]] = likelihoodOfSNV(value, depth, pUsed1, overdispersionValue)
      pListChanged=T
    }
    firstLikelihood = (pList[[as.character(pUsed1)]])
    if (!as.character(pUsed2) %in% names(pList)) {
      pList[[as.character(pUsed2)]] = (likelihoodOfSNV(value, depth, pUsed2, overdispersionValue))
      pListChanged=T
    }
    secondLikelihood = (pList[[as.character(pUsed2)]])
    finalLikelihood = log(0.5 * firstLikelihood + 0.5 * secondLikelihood)
  } else if (cn == 0) {
    finalLikelihood = log(likelihoodOfSNV(value, depth, multiplierDueToMapping * 0.5, overdispersionValue))
  } else if (stateUsed %in% c("LOH", "LOHDup")) {
    if (cn == 2) {
     pUsed1=round(multiplierDueToMapping * max(0.01, 0.5 - pur/2), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, 0.5 + pur/2), digits=2)
    } else if (cn == 3) {
      pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 2 * pur) / (2 + pur)), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, (1 - pur) / (2 + pur)), digits=2)
    } else if (cn == 4) {
      pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 3 * pur) / (2 + 2 * pur)), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, (1 - pur) / (2 + 2 * pur)), digits=2)
    }
    if (!as.character(pUsed1) %in% names(pList)) {
      pList[[as.character(pUsed1)]] = likelihoodOfSNV(value, depth, pUsed1, overdispersionValue)
      pListChanged=T
    }
    if (!as.character(pUsed2) %in% names(pList)) {
      pList[[as.character(pUsed2)]] = (likelihoodOfSNV(value, depth, pUsed2, overdispersionValue))
      pListChanged=T
    }
    finalLikelihood = log(0.5 * pList[[as.character(pUsed1)]] + 
                            0.5 * pList[[as.character(pUsed2)]])
  }  else if (stateUsed == "CNVboth" | stateUsed == "normal") {
    pUsed = round((multiplierDueToMapping * 0.5), digits=2)
    if (!as.character(pUsed) %in% names(pList)) {
      pList[[as.character(pUsed)]] = (likelihoodOfSNV(value, depth, pUsed, overdispersionValue))
      pListChanged=T
    }
    finalLikelihood = log(pList[[as.character(pUsed)]])
  } else if (stateUsed == "CNVcomplex2" | stateUsed == "CNVcomplex3") {
    if (cn == 5) {
      pUsed1=round(multiplierDueToMapping * max(0.01, (1 + pur) / (2 + 3 * pur)), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 2 * pur) / (2 + 3* pur)), digits=2)

    } else if (cn == 6) {
      pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 1 * pur) / (2 + 4 * pur)), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 3 * pur) / (2 + 4 * pur)), digits=2)
    } else if (cn == 7) {
      if (stateUsed == "CNVcomplex3") {
        pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 3 * pur) / (2 + 5 * pur)), digits=2)
        pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 2 * pur) / (2 + 5 * pur)), digits=2)
      } else {
        pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 4 * pur) / (2 + 5 * pur)), digits=2)
        pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 1 * pur) / (2 + 5 * pur)), digits=2)
    }
    } else if (cn == 8) {
      if (stateUsed == "CNVcomplex3") {
        pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 4 * pur) / (2 + 6 * pur)), digits=2)
        pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 2 * pur) / (2 + 6 * pur)), digits=2)
      } else {
        pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 1 * pur) / (2 + 6 * pur)), digits=2)
        pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 5 * pur) / (2 + 6 * pur)), digits=2)
      }
    }
    if (!as.character(pUsed1) %in% names(pList)) {
      pList[[as.character(pUsed1)]] = likelihoodOfSNV(value, depth, pUsed1, overdispersionValue)
      pListChanged=T
    }
    if (!as.character(pUsed2) %in% names(pList)) {
      pList[[as.character(pUsed2)]] = (likelihoodOfSNV(value, depth, pUsed2, overdispersionValue))
      pListChanged=T
    }
    finalLikelihood = log(0.5 * pList[[as.character(pUsed1)]] + 
                            0.5 * pList[[as.character(pUsed2)]])
  } else {
    print("State used was not described by previous condition")
    print(pur)
    print(cn)
    print(stateUsed)
    stop()
  }
  if (pListChanged){
    return(list(finalLikelihood, pListChanged, pList))} else {
      return(list(finalLikelihood, pListChanged))
    }
}


extractVariancesFromBAF <- function(bafTable, expectedValue) {
  freqs = as.numeric(bafTable[,5])
  depths = as.numeric(bafTable[,6])
  depths = depths[which(freqs > expectedValue - 0.1 & freqs < expectedValue + 0.1)]
  freqs = freqs[which(freqs > expectedValue - 0.1 & freqs < expectedValue + 0.1)]
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
    rollingLength = max(10, round(length(forEstimationFreqs) / 10))
    overdispersion = rep(0, length(forEstimationFreqs))
    for (i in 1:length(vecOfSDs)) {
      vecOfSDs[i] = Qn(forEstimationFreqs[max(1,i-rollingLength):min(length(vecOfSDs), i+rollingLength)] * forEstimationDepths[i]) ** 2
      varPredictedByBinom = forEstimationDepths[i] * expectedValue * (1 - expectedValue)
      overdispersion[i] = (vecOfSDs[i]) / varPredictedByBinom
    }
    overdispersion[which(overdispersion < 1)] = 1
    overdispersion[which(overdispersion >4)] = 4
    #plot(overdispersion ~ forEstimationDepths)
    return(cbind(overdispersion, forEstimationDepths))
  }
}
