
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

determineHeterozygousPositions <- function(freq, depth) {
  prob = pbinom(round(freq * depth), prob=0.48, size=depth)
  if (prob > 0.01 & prob < 0.99) {
    return(T)
  } else {
    return(F)
  }
}



likelihoodOfSNVBasedOnCN <- function(value, depth, pur, cn) {
  multiplierDueToMapping = 0.475 / 0.5
  overallNumberOfReads <- (1 - pur) * 2 + pur * (cn)
  if (cn != 0 & cn != 2) {
    numberOfReadsSupportiveOne1 <- (1 - pur) + pur * (cn - 1)
    numberOfReadsSupportiveOne2 <- 1 # (1 - pur) + pur * 1
    finalLikelihood = log(0.5 * likelihoodOfSNV(value, depth, multiplierDueToMapping * min(0.999, numberOfReadsSupportiveOne1 / overallNumberOfReads)) + 
                       0.5 * likelihoodOfSNV(value, depth, multiplierDueToMapping * max(0.001, numberOfReadsSupportiveOne2 / overallNumberOfReads))
                       )
  } else if (cn == 0) {
    finalLikelihood = log(likelihoodOfSNV(value, depth, 0.5))
  } else if (cn == 2) {
    if (pur == 0) {
      finalLikelihood = log(likelihoodOfSNV(value, depth, multiplierDueToMapping * 0.5))
    } else {
      finalLikelihood = log(0.5 * likelihoodOfSNV(value, depth, multiplierDueToMapping * max(0.001, 0.5 - pur/2)) + 
                              0.5 * likelihoodOfSNV(value, depth, multiplierDueToMapping * min(0.999, 0.5 + pur/2)))
    }
  }
  return(finalLikelihood)
}
