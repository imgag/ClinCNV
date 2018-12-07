
likelihoodOfSNV <- function(a, b, p) {

    value = dbinom(a, size=b, prob=p)
    
  if (is.nan(value) | value < 10**-30) return(10**-30)
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

determineHeterozygousPositions <- function(freq, depth, probAB=0.48) {
  prob = pbinom(round(freq * depth), prob=probAB, size=depth)
  if ((prob > 0.01 & prob < 0.99)) {
    return(T)
  } else {
    return(F)
  }
}



likelihoodOfSNVBasedOnCN <- function(value, depth, pur, cn, stateUsed, multiplierOfSNVsDueToMapping, pList) {
  multiplierDueToMapping = multiplierOfSNVsDueToMapping / 0.5
  overallNumberOfReads <- (1 - pur) * 2 + pur * (cn)
  pListChanged = F
  if (cn != 0 & cn != 2 & stateUsed == "CNV") {
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
      pList[[as.character(pUsed1)]] = likelihoodOfSNV(value, depth, pUsed1)
      pListChanged=T
    }
    if (!as.character(pUsed2) %in% names(pList)) {
      pList[[as.character(pUsed2)]] = (likelihoodOfSNV(value, depth, pUsed2))
      pListChanged=T
    }
    finalLikelihood = log(0.5 * pList[[as.character(pUsed1)]] + 
                            0.5 * pList[[as.character(pUsed2)]])
  }  else if (stateUsed == "CNVboth" | stateUsed == "normal") {
    pUsed = round((multiplierDueToMapping * 0.5), digits=2)
    if (!as.character(pUsed) %in% names(pList)) {
      pList[[as.character(pUsed)]] = (likelihoodOfSNV(value, depth, pUsed))
      pListChanged=T
    }
    finalLikelihood = log(pList[[as.character(pUsed)]])
  } else if (stateUsed == "CNVcomplex") {
    if (cn == 5) {
      pUsed1=round(multiplierDueToMapping * max(0.01, (1 + pur) / (2 + 3 * pur)), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 2 * pur) / (2 + 3* pur)), digits=2)
      pUsed3=pUsed1
      pUsed4=pUsed2
    } else if (cn == 6) {
      pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 1 * pur) / (2 + 4 * pur)), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 3 * pur) / (2 + 4 * pur)), digits=2)
      pUsed3=pUsed1
      pUsed4=pUsed2
    } else if (cn == 7) {
      pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 3 * pur) / (2 + 5 * pur)), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 2 * pur) / (2 + 5 * pur)), digits=2)
      pUsed3=round(multiplierDueToMapping * max(0.01, (1 + 4 * pur) / (2 + 5 * pur)), digits=2)
      pUsed4=round(multiplierDueToMapping * min(0.99, (1 + 1 * pur) / (2 + 5 * pur)), digits=2)
    } else if (cn == 8) {
      pUsed1=round(multiplierDueToMapping * max(0.01, (1 + 4 * pur) / (2 + 6 * pur)), digits=2)
      pUsed2=round(multiplierDueToMapping * min(0.99, (1 + 2 * pur) / (2 + 6 * pur)), digits=2)
      pUsed3=round(multiplierDueToMapping * max(0.01, (1 + 1 * pur) / (2 + 6 * pur)), digits=2)
      pUsed4=round(multiplierDueToMapping * min(0.99, (1 + 5 * pur) / (2 + 6 * pur)), digits=2)
    }
    if (!as.character(pUsed1) %in% names(pList)) {
      pList[[as.character(pUsed1)]] = likelihoodOfSNV(value, depth, pUsed1)
      pListChanged=T
    }
    if (!as.character(pUsed2) %in% names(pList)) {
      pList[[as.character(pUsed2)]] = (likelihoodOfSNV(value, depth, pUsed2))
      pListChanged=T
    }
    if (!as.character(pUsed3) %in% names(pList)) {
      pList[[as.character(pUsed3)]] = likelihoodOfSNV(value, depth, pUsed3)
      pListChanged=T
    }
    if (!as.character(pUsed4) %in% names(pList)) {
      pList[[as.character(pUsed4)]] = (likelihoodOfSNV(value, depth, pUsed4))
      pListChanged=T
    }
    finalLikelihood = log(0.25 * pList[[as.character(pUsed1)]] + 
                            0.25 * pList[[as.character(pUsed2)]] +
                            0.25 * pList[[as.character(pUsed3)]] + 
                            0.25 * pList[[as.character(pUsed4)]])
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
