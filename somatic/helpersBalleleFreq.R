
likelihoodOfSNV <- function(a, b, p) {
  if ((b * p >= 10 & b * (1 - p) >= 10) ) {
    value = return_norm_likelik((a - b * p) / sqrt(b * p * (1-p)))
    } else {
    value = dbinom(a, size=b, prob=p)
    }
  if (is.nan(value) | value < 10**-40) return(10**-40)
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

determineHeterozygousPositions <- function(freq, depth) {
  prob = pbinom(round(freq * depth), prob=0.48, size=depth)
  if (prob > 0.01 & prob < 0.99) {
    return(T)
  } else {
    return(F)
  }
}



likelihoodOfSNVBasedOnCN <- function(value, depth, pur, cn, pList) {
  multiplierDueToMapping = 0.48 / 0.5
  overallNumberOfReads <- (1 - pur) * 2 + pur * (cn)
  pListChanged = F
  if (cn != 0 & cn != 2) {
    numberOfReadsSupportiveOne1 <- (1 - pur) + pur * (cn - 1)
    numberOfReadsSupportiveOne2 <- 1 # (1 - pur) + pur * 1
    pUsed1=round(multiplierDueToMapping * min(0.99, numberOfReadsSupportiveOne1 / overallNumberOfReads), digits=2)
    pUsed2=round(multiplierDueToMapping * max(0.01, numberOfReadsSupportiveOne2 / overallNumberOfReads), digits=2)
    if (!as.character(pUsed1) %in% names(pList)) {
      pList[[as.character(pUsed1)]] = likelihoodOfSNV(value, depth, pUsed1)
      pListChanged=T
    }
    firstLikelihood = (pList[[as.character(pUsed1)]])
    if (!as.character(pUsed2) %in% names(pList)) {
      pList[[as.character(pUsed2)]] = (likelihoodOfSNV(value, depth, pUsed2))
      pListChanged=T
    }
    secondLikelihood = (pList[[as.character(pUsed2)]])
    finalLikelihood = log(0.5 * firstLikelihood + 0.5 * secondLikelihood)
  } else if (cn == 0) {
    finalLikelihood = log(likelihoodOfSNV(value, depth, multiplierDueToMapping * 0.5))
  } else if (cn == 2) {
    if (pur == 0) {
      pUsed = round((multiplierDueToMapping * 0.5), digits=2)
      if (!as.character(pUsed) %in% names(pList)) {
        pList[[as.character(pUsed)]] = (likelihoodOfSNV(value, depth, pUsed))
        pListChanged=T
      }
      finalLikelihood = log(pList[[as.character(pUsed)]])
    } else {
      pUsed1=round(multiplierDueToMapping * max(0.01, 0.5 - pur/2), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, 0.5 + pur/2), digits=2)
      if (!as.character(pUsed1) %in% names(pList)) {
        pList[[as.character(pUsed1)]] = likelihoodOfSNV(value, depth, pUsed1)
        pListChanged=T
      }
      firstLikelihood = (pList[[as.character(pUsed1)]])
      if (!as.character(pUsed2) %in% names(pList)) {
        pList[[as.character(pUsed2)]] = (likelihoodOfSNV(value, depth, pUsed2))
        pListChanged=T
      }
      finalLikelihood = log(0.5 * pList[[as.character(pUsed1)]] + 
                              0.5 * pList[[as.character(pUsed2)]])
    }
  }
  if (pListChanged){
  return(list(finalLikelihood, pListChanged, pList))} else {
    return(list(finalLikelihood, pListChanged))
  }
}
