
EstimateModeSimple <- function(x) {

  density_of_x <-  density(x, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  mu
}



likelihoodOfGaussianMixture <- function(location, value , sd_to_start, robustPercentage, minSizeOfCluster, esimtatedVarianceFromSampleNoise, lowerBoundSD) {
  initLocation = location
  if (length(location) == 1) {
    cluster_weights = rep(1/ length(location), length(location))
    points_likeliks <- log(sapply(1:length(location), 
                                  function(i) {
                                    data_for_mean = (value - location[i]) / (sd_to_start * esimtatedVarianceFromSampleNoise); 
                                    return(cluster_weights[i] * 
                                             with(as.data.frame(data_for_mean), 
                                                  return_likelik(data_for_mean)) / (sd_to_start * esimtatedVarianceFromSampleNoise))
                                  }))
    minLikelik = quantile(points_likeliks, robustPercentage)
    points_likeliks[which(points_likeliks < minLikelik)] = minLikelik
    return(sum(points_likeliks))
  } else {
    sd_counter = sd_to_start
    eps = 0.01
    cluster_weights = rep(1/ length(location), length(location))
    points_likeliks <- (sapply(1:length(location), 
                               function(i) {
                                 data_for_mean = (value - location[i]) / (sd_counter * esimtatedVarianceFromSampleNoise); 
                                 return(cluster_weights[i] * 
                                          with(as.data.frame(data_for_mean), 
                                               return_likelik(data_for_mean)) / (sd_counter * esimtatedVarianceFromSampleNoise))
                               }))
    robust_res <- list(rowSums(points_likeliks))
    vect_sum <- robust_res[[1]]
    previous_loglik <- sum(log(vect_sum))
    # FIRST LIKELIHOOD FIXED, NOW WE MOVE TO ITERATIONS
    delta_loglik <- eps * 2
    counter=0
    while (delta_loglik > eps) {
      weights = points_likeliks * (1 / vect_sum) / length(value)
      cluster_weights <- colSums(weights) + 10**-100
      locationTmp = location
      for (l in 1:length(locationTmp)) {
        locationTmp[l] = sum(weights[,l] * value) / cluster_weights[l]
      }
      multiplierForLocation = locationTmp[which.max(cluster_weights[l])] / location[which.max(cluster_weights[l])]
      if (abs(log2(initLocation / (location * multiplierForLocation))) < 0.1) {
        location = location * multiplierForLocation
      } 
      sd_counter_tmp = sqrt(sum(sapply(1:length(location), 
                                       function(i) {
                                         weights[,i] * (value - location[i])^2
                                       })))
      if (sd_counter_tmp > lowerBoundSD) {
        sd_counter = sd_counter_tmp
      }
      #sd_to_start = getSdsAllSamples(sd_counter)
      points_likeliks <- (sapply(1:length(location), 
                                 function(i) {
                                   data_for_mean = (value - location[i]) / (sd_counter * esimtatedVarianceFromSampleNoise); 
                                   return(cluster_weights[i] * 
                                            with(as.data.frame(data_for_mean), 
                                                 return_likelik(data_for_mean)) / (sd_counter * esimtatedVarianceFromSampleNoise))
                                 }))
      robust_res <- list(rowSums(points_likeliks))
      vect_sum <- robust_res[[1]]
      current_loglik <- sum(log(vect_sum))
      counter = counter + 1
      delta_loglik = current_loglik - previous_loglik
      previous_loglik = current_loglik
    }
    if (counter > 500)    print(paste("We finished fitting in", counter, "steps"))
    weights = points_likeliks * (1 / vect_sum) / length(value)
    cluster_weights <- colSums(weights)
    # make small values as nulls
    cluster_weights[which(cluster_weights < minSizeOfCluster / length(value))] = 0
    points_likeliks <- (sapply(1:length(location), 
                               function(i) {
                                 data_for_mean = (value - location[i]) / (sd_counter * esimtatedVarianceFromSampleNoise); 
                                 return(cluster_weights[i] * 
                                          with(as.data.frame(data_for_mean), 
                                               return_likelik(data_for_mean)) / (sd_counter * esimtatedVarianceFromSampleNoise))
                               }))
    robust_res <- list(rowSums(points_likeliks))
    vect_sum <- robust_res[[1]]
    current_loglik <- sum(log(vect_sum))
    weights = points_likeliks * (1 / vect_sum) / length(value)
    cluster_weights <- colSums(weights)
    
    
  }
  
  return(list(current_loglik, cluster_weights / sum(cluster_weights), sd_counter))
}

checkConnectivity = function(covOne, covTwo) {
  whichBothNonHomo = which(covOne > 0.5 & covTwo > 0.5)
  newlm1 = rlm(covOne[whichBothNonHomo], covTwo[whichBothNonHomo])
  newlm2 = rlm(covTwo[whichBothNonHomo], covOne[whichBothNonHomo])
  QnResid1 = Qn(newlm1$residuals)
  QnResid2 = Qn(newlm2$residuals)
  if (length(which(newlm1$residuals > 2.5 * QnResid1)) < 0.05 * length(whichBothNonHomo) | length(which(newlm2$residuals > 2.5 * QnResid2)) < 0.05 * length(whichBothNonHomo)) {
    return(T)
  }
  plot(covOne ~ covTwo, main="CHECKED")
  return(F)
}

checkConnectivityMed = function(covOne, covTwo, sampleVariability) {
  whichBothNonHomo = which(covOne > 0.3 & covTwo > 0.3)
  if (length(whichBothNonHomo) < 5) {
    return(F)
  }
  angle = median(covOne / covTwo)
  distances <- sapply(1:length(whichBothNonHomo), function(i) {y0  = covOne[whichBothNonHomo[i]]; 
  x0 = covTwo[whichBothNonHomo[i]];
  (-angle * x0 + y0) / sqrt(angle ** 2 + 1)
  })
  if (length(which(abs(distances) < sqrt(3/2) - 1)) < 5) {
    return(F)
  }
  QnDist = Qn(distances[which(abs(distances) < sqrt(3/2) - 1)]) * sampleVariability
  if (length(which(distances > 3 * QnDist)) < 0.05 * length(whichBothNonHomo)) {
    return(T)
  }
  plot(covOne ~ covTwo, main="CHECKED")
  return(F)
}

checkConnectivityComplex = function(j) {
  if (j == length(localSdsOfProbes)) {
    return(F)
  }
  sdFirst=localSdsOfProbes[j]
  sdSecond = localSdsOfProbes[j+1]
  sdThree = NULL
  if (j + 2 <= length(localSdsOfProbes)) {
    sdThree = localSdsOfProbes[j+2]
  }
  locations1 = bestLocations[[j]]
  locations2 = bestLocations[[j+1]]
  if (j + 2 <= length(localSdsOfProbes)) {
    locations3 = bestLocations[[j+2]]
  }
  weightsOne = bestWeights[[j]]
  weightsTwo = bestWeights[[j+1]]
  if (j + 2 <= length(localSdsOfProbes)) {
    weightsThree = bestWeights[[j+2]]
  }
  if (is.null(weightsOne) | is.null(weightsTwo)) {
    return(F)
  }
  sdsMultipliers = multipliersSamples
  dataOne = coverageToWorkWith[j,]
  dataOne[which(dataOne < min(locations1))] = min(locations1) + rnorm(length(which(dataOne < min(locations1))), sd=sdFirst)
  dataOne[which(dataOne > max(locations1))] = max(locations1) + rnorm(length(which(dataOne > max(locations1))), sd=sdFirst)
  dataTwo = coverageToWorkWith[j+1,]
  dataTwo[which(dataTwo < min(locations2))] = min(locations2) + rnorm(length(which(dataTwo < min(locations2))), sd=sdSecond)
  dataTwo[which(dataTwo > max(locations2))] = max(locations2) + rnorm(length(which(dataTwo > max(locations2))), sd=sdSecond)
  if (j + 2 <= length(localSdsOfProbes)) {
    dataThree = coverageToWorkWith[j+2,]
    dataThree[which(dataThree < min(locations3))] = min(locations3) + rnorm(length(which(dataThree < min(locations3))), sd=sdThree)
    dataThree[which(dataThree > max(locations3))] = max(locations3) + rnorm(length(which(dataThree > max(locations3))), sd=sdThree)
  }
  minSizeOfCluster = 2
  resultList = likelihoodForTwoNeighbors(sdFirst, sdSecond, locations1, locations2, weightsOne, weightsTwo, sdsMultipliers, dataOne, dataTwo, minSizeOfCluster)
  if (length(resultList) == 1) {
    return(F)
  }
  clusterWeights = resultList[[2]]
  arrayOfDiags <- c(sum(diag(clusterWeights)))
  for (i in 1:2) {
    arrayOfDiags = c(arrayOfDiags, sum(clusterWeights[,-i]))
    arrayOfDiags = c(arrayOfDiags, sum(clusterWeights[-i,]))
  }
  if (max(arrayOfDiags) > 0.95) {
    return(T)
  } else {
    if (!is.null(weightsThree)) {
      resultListAdd = likelihoodForTwoNeighbors(sdFirst, sdThree, locations1, locations3, weightsOne, weightsThree, sdsMultipliers, dataOne, dataThree, minSizeOfCluster)
    }
    if (!is.null(weightsThree)) {
      clusterWeightsNew = resultListAdd[[2]]
      arrayOfDiagsNew <- c(sum(diag(clusterWeightsNew)))
      for (i in 1:2) {
        arrayOfDiagsNew = c(arrayOfDiagsNew, sum(clusterWeightsNew[,-i]))
        arrayOfDiagsNew = c(arrayOfDiagsNew, sum(clusterWeightsNew[-i,]))
      }
      if (max(arrayOfDiagsNew) > 0.95) {
        return(T)
      } 
    }
    return(F)
  }
}


findFinalState <- function(coverageNeededToCheck, toyBedFilePolymorphCurrent, multipliersSamples, numberOfClusterAnalysed, plotting, chromX = F, folder_name_mcnv) {
  startOfmCNV = 2
  endOfmCNV = (nrow(coverageNeededToCheck) - 1)
  if (toyBedFilePolymorphCurrent[1,3] - toyBedFilePolymorphCurrent[1,2] < 250 | nrow(coverageNeededToCheck) < 3) {
    startOfmCNV = 1
  }
  if (toyBedFilePolymorphCurrent[nrow(coverageNeededToCheck),3] - toyBedFilePolymorphCurrent[nrow(coverageNeededToCheck),2] < 250 | nrow(coverageNeededToCheck) < 3) {
    endOfmCNV = nrow(coverageNeededToCheck)
  }
  #if (nrow(coverageNeededToCheck) > 2) {
    coverageSummarised = apply(coverageNeededToCheck[startOfmCNV:endOfmCNV,,drop=F], 2, median)
  #}
    coverageSummarised = coverageSummarised / quantile(coverageSummarised, 0.8)
  notHomozygousDeletions = which(coverageSummarised >= 0.5)
  if (length(which(coverageSummarised <= 0.25)) > 0) {
    homozygousDelShit = median(coverageSummarised[which(coverageSummarised <= 0.25)]) ** 2
    coverageSummarised = sqrt(abs(coverageSummarised ** 2 - homozygousDelShit))
  }
  for (i in 1:20) {
    locations[[i]] = sqrt(1:30/i)
  }
  modeOfCovSummarised = EstimateModeSimple(coverageSummarised[notHomozygousDeletions])
  coverageSummarised = coverageSummarised / modeOfCovSummarised
  notHomozygousDeletions = which(coverageSummarised >= 0.45)
  bestLoc = sqrt(0:20/2)
  sdNormalised = Qn(coverageSummarised[which(coverageSummarised > sqrt(1/2) & coverageSummarised < sqrt(3/2))])
  #likelikAndWeights = likelihoodOfGaussianMixture(1, coverageSummarised[notHomozygousDeletions], sdNormalised, 
  #                                                0.05 * (length(notHomozygousDeletions)) / length(coverageSummarised), 0.1, 
  #                                                multipliersSamples[notHomozygousDeletions], lowerBoundOfSD)
  bestLikelik = 10**100
  bestSD = NULL
  bestDivisor = 2
  bestWeight = c(1)
  bestLoc = c(1)
  #possibleLocations = round(medianOfCoverage ** 2 * 2)
  #possibleLocations = unique(c(possibleLocations, possibleLocations + -5:5))
  coverageSummarisedCleaned = coverageSummarised[which(multipliersSamples < quantile(multipliersSamples, 0.95))]
  notHomozygousDeletionsCleaned = which(coverageSummarisedCleaned >= 0.45)
  for (j in 1:length(locations)) {
    #if (!j %in% possibleLocations) next
    vecOfMeans = locations[[j]]
    vecOfMeans = vecOfMeans[which(vecOfMeans > min(coverageSummarisedCleaned[notHomozygousDeletionsCleaned]) - 0.25 & vecOfMeans < max(coverageSummarisedCleaned[notHomozygousDeletionsCleaned]) + 0.25)]
    if (length(which(coverageSummarisedCleaned <= 0.25)) / length(coverageSummarisedCleaned) > 0.01) {
      vecOfMeans = unique(round(c(vecOfMeans, 1 / sqrt(j)) * 100000))
      vecOfMeans = vecOfMeans / 100000
    }
    if (length(vecOfMeans) == 1) next
    likelikAndWeights = likelihoodOfGaussianMixture(vecOfMeans, coverageSummarisedCleaned[notHomozygousDeletionsCleaned], sdNormalised, 
                                                    0.05 * (length(notHomozygousDeletionsCleaned)) / length(coverageSummarisedCleaned), 0.1, 
                                                    multipliersSamples[notHomozygousDeletionsCleaned], lowerBoundOfSD)
    firstSignifCluster = min(which(likelikAndWeights[[2]] > 0.025))
    if (length(notHomozygousDeletionsCleaned) < length(coverageSummarisedCleaned)) {
      firstSignifCluster = 1
    }
    lastSignifCluster = max(which(likelikAndWeights[[2]] > 0.025))
    if (firstSignifCluster < lastSignifCluster - 1)
      if (min(likelikAndWeights[[2]][(firstSignifCluster + 1):(lastSignifCluster - 1)]) < 0.01) next
    if (length(which(coverageSummarisedCleaned <= 0.25)) / length(coverageSummarisedCleaned) > 0.01) {
      whichHomoDel = which.min(vecOfMeans - 1/sqrt(j))
      if (likelikAndWeights[[2]][whichHomoDel] < 0.01) {
        next
      }
    }
    tmpLikelik =  (
      -2 * likelikAndWeights[[1]] + (lastSignifCluster - firstSignifCluster  + 1) * log(length(notHomozygousDeletionsCleaned))
    )

    bestWeightCheckForEvenDominance = c(likelikAndWeights[[2]] * length(notHomozygousDeletionsCleaned)) + 0.1
    bestWeightCheckForEvenDominance = bestWeightCheckForEvenDominance / sum(bestWeightCheckForEvenDominance)
    bestLocEvenDominance = c(vecOfMeans)
    bestLocEvenDominance = round(bestLocEvenDominance ** 2 * j)
    if (!chromX & (sum(bestWeightCheckForEvenDominance[which(bestLocEvenDominance %% 2 == 0)]) * length(notHomozygousDeletionsCleaned) + length(which(coverageSummarisedCleaned < 0.5))) / length(coverageSummarisedCleaned) < 0.4) next

    if (tmpLikelik < bestLikelik | is.null(bestSD)) {
      if (3 * likelikAndWeights[[3]] < 1 - sqrt((j-1)/j)) {
        bestLikelik= tmpLikelik
        bestLoc = vecOfMeans
        bestSD = likelikAndWeights[[3]]
        bestWeight = likelikAndWeights[[2]]
        bestDivisor = j
      }
    }
  }
  if (length(bestWeight) == 1) {
    bestLoc = sqrt(1:20/2)
    bestWeight = rep(0, 20)
    for (c in 1:length(coverageSummarised)) {
      if (coverageSummarised[c] > 0.5) {
        closestLoc = which.min(abs(coverageSummarised[c] - bestLoc))
        bestWeight[closestLoc] = bestWeight[closestLoc] + 1
      }
    }
    bestWeight = bestWeight / sum(bestWeight)
    bestSD = sdNormalised
  }
  if (is.null(bestSD)) {
    bestSD = sdNormalised
  }
  bestLoc = unique(c(0, bestLoc))
  if (length(which(coverageSummarised <= 0.25)) > 0) {
    bestLoc[which(bestLoc ==0)] = median(coverageSummarised[which(coverageSummarised <= 0.25)])
  }
  bestWeight = c(length(coverageSummarised) - length(notHomozygousDeletions), bestWeight * length(notHomozygousDeletions)) + 0.1
  bestWeight = bestWeight / sum(bestWeight)
  
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  copy_number_likeliks <- abs(sweep(rep.row(coverageSummarised, length(bestLoc)), 1, bestLoc, FUN="-"))
  copy_number_likeliks <- t(apply(copy_number_likeliks, 1, function(x){dnorm(x, sd = bestSD)}))
  copy_number_likeliks <- sweep(copy_number_likeliks, 1, bestWeight, FUN="*")
  copy_number <- bestLoc[apply(copy_number_likeliks, 2, which.max)]
  copy_number = round(copy_number ** 2 * bestDivisor)
  coloursP = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414, 337)]
  coloursP = c(coloursP, coloursP)
  diagnosticPlot = (length(which(copy_number != as.numeric(names(sort(table(copy_number),decreasing=TRUE)[1])))) >= 0.05 * length(copy_number))
  if (diagnosticPlot == T & plotting) {
    fileName = paste(numberOfClusterAnalysed, toyBedFilePolymorphCurrent[1,1], toyBedFilePolymorphCurrent[1,2], toyBedFilePolymorphCurrent[nrow(toyBedFilePolymorphCurrent),3], sep="_")
    png(paste0(folder_name_mcnv, fileName, ".png"), width=length(copy_number) * 3, height=800)
    plot(coverageSummarised ** 2 * bestDivisor, col="black", pch=21,bg=coloursP[(copy_number + 1)])
    abline(h=bestLoc ** 2 * bestDivisor)
    dev.off()
  }
  #plot(density(coverageSummarised, bw='SJ'))
  return(copy_number)
}





likelihoodForTwoNeighbors = function(sdFirst, sdSecond, locations1, locations2, weightsOne, weightsTwo, sdsMultipliers, dataOne, dataTwo, minSizeOfCluster) {
  covarianceMatrix = matrix(c(sdFirst**2, 0, 0, sdSecond**2 ), nrow=2)
  matrixOfLikeliks = array(0, c(length(dataOne), length(locations1), length(locations2)))
  for (k in 1:length(dataOne)) {
    for (i in 1:length(locations1)) {
      for (j in 1:length(locations2)) {
        matrixOfLikeliks[k, i, j] = exp( -0.5 * 
                                           (
                                             log( covarianceMatrix[1,1] * (sdsMultipliers[k] ** 4) * covarianceMatrix[2,2] ) + 
                                               ((dataOne[k] - locations1[i]) ** 2 * 1 / ((sdsMultipliers[k]) ** 2 * covarianceMatrix[1,1])) + 
                                               ((dataTwo[k] - locations2[j]) ** 2 * 1 / ((sdsMultipliers[k]) ** 2 * covarianceMatrix[2,2])) +
                                               2 * log( 2 * pi)
                                           )
        )
        if (is.na(matrixOfLikeliks[k, i, j])) matrixOfLikeliks[k, i, j] = 0
      }
    }
  }
  
  
  clusterWeights <- diag(1 / length(locations1), nrow=length(locations1), ncol=length(locations2))
  for (i in 1:length(locations1)) {
    for (j in 1:length(locations2)) {
      clusterWeights[i,j] = weightsOne[i] * weightsTwo[j]
    }
  }
  clusterWeights = clusterWeights + min(10**-10, 1/((length(locations1) ** 2)))
  clusterWeights = clusterWeights / sum(clusterWeights)
  
  
  gammaNK = array(0, c(length(dataOne), length(locations1), length(locations2)))
  for (k in 1:length(dataOne)) {
    gammaNK[k,,] = clusterWeights * matrixOfLikeliks[k,,]
  }
  
  eachPointLoglik = rep(0, length(dataOne))
  for (k in 1:length(dataOne)) {
    eachPointLoglik[k] = log(sum(clusterWeights * matrixOfLikeliks[k,,]))
    if (is.na(eachPointLoglik[k]) | is.infinite(eachPointLoglik[k])) eachPointLoglik[k] = -10**100
  }
  previousLoglik = sum(eachPointLoglik)
  
  epsilon = 0.1
  deltaLoglik = 2 * epsilon
  counter = 1
  while (epsilon < deltaLoglik) {
    print(paste("Iteration number", counter))
    counter = counter + 1
    print(paste("Change in likelihood", deltaLoglik))
    for (k in 1:length(dataOne)) {
      clusterWeightIJ = sum(gammaNK[k,,])
      if (is.nan(clusterWeightIJ)) {
        return(F)
        stop("Numerical error of EM algorithm")
      }
      if (!is.na(clusterWeightIJ) | clusterWeightIJ > 0) {
        gammaNK[k,,] = gammaNK[k,,] / (clusterWeightIJ)
      } else {
        print("THIS IS THERO! FUCK")
      }
    }
    for (i in 1:length(locations1)) {
      for (j in 1:length(locations2)) {
        clusterWeights[i,j] = sum(gammaNK[,i,j])
      }
    }
    clusterWeights = clusterWeights / length(dataOne)
    
    
    eachPointLoglik = rep(0, length(dataOne))
    for (k in 1:length(dataOne)) {
      eachPointLoglik[k] = log(sum(clusterWeights * matrixOfLikeliks[k,,]))
      if (is.na(eachPointLoglik[k])) eachPointLoglik[k] = -10**100
    }
    currentLoglik = sum(eachPointLoglik)
    deltaLoglik = currentLoglik - previousLoglik
    previousLoglik = currentLoglik
    if (deltaLoglik > epsilon) {
      gammaNK = array(0, c(length(dataOne), length(locations1), length(locations2)))
      for (k in 1:length(dataOne)) {
        gammaNK[k, , ] = clusterWeights * matrixOfLikeliks[k, , ]
      }
    } else {
      #return(list(clusterWeights, currentLoglik))
      clustersToExclude <- which(clusterWeights < minSizeOfCluster / length(dataOne))
      clusterWeights[clustersToExclude] = 0
      gammaNK = array(0, c(length(dataOne), length(locations1), length(locations2)))
      for (k in 1:length(dataOne)) {
        gammaNK[k, , ] = clusterWeights * matrixOfLikeliks[k, , ]
      }
      for (k in 1:length(dataOne)) {
        clusterWeightIJ = sum(gammaNK[k,,]) 
        if (!is.na(clusterWeightIJ) | clusterWeightIJ > 0) {
          gammaNK[k,,] = gammaNK[k,,] / (clusterWeightIJ)
        } else {
          print("THIS IS THERO! FUCK")
        }
      }
      for (i in 1:length(locations1)) {
        for (j in 1:length(locations2)) {
          clusterWeights[i,j] = sum(gammaNK[,i,j])
        }
      }
      clusterWeights = clusterWeights / length(dataOne)
      eachPointLoglik = rep(0, length(dataOne))
      for (k in 1:length(dataOne)) {
        eachPointLoglik[k] = log(sum(clusterWeights * matrixOfLikeliks[k,,]))
        if (is.na(eachPointLoglik[k])) eachPointLoglik[k] = -10**100
      }
      currentLoglik = sum(eachPointLoglik)
      #return(list(clusterWeights, currentLoglik))
      
      oldBIC = -2 * currentLoglik + log(length(dataOne)) * length(which(clusterWeights != 0))
      #print(oldBIC)
      while (length(which(clusterWeights != 0)) > 0) {
        #print(paste("BIC iteration, BIC is equal to", oldBIC))
        print(paste("Number of non zero clusters", length(which(clusterWeights > 0))))
        if (length(which(clusterWeights > 0)) <= 1){
          print(length(which(clusterWeights > 0)) <= 1)
          return(list(100, matrix(c(0,0,0,1),nrow=2)))
        }
        
        tmpClusterWeights = clusterWeights
        minClusterWeight = min(tmpClusterWeights[which(tmpClusterWeights > 0)])
        tmpClusterWeights[which(tmpClusterWeights <= minClusterWeight)] = 0
        gammaNK = array(0, c(length(dataOne), length(locations1), length(locations2)))
        for (k in 1:length(dataOne)) {
          gammaNK[k, , ] = tmpClusterWeights * matrixOfLikeliks[k, , ]
        }
        for (k in 1:length(dataOne)) {
          clusterWeightIJ = sum(gammaNK[k,,]) 
          if (!is.na(clusterWeightIJ) | clusterWeightIJ > 0) {
            gammaNK[k,,] = gammaNK[k,,] / (clusterWeightIJ)
          } else {
            print("THIS IS THERO! FUCK")
          }
        }
        for (i in 1:length(locations1)) {
          for (j in 1:length(locations2)) {
            tmpClusterWeights[i,j] = sum(gammaNK[,i,j])
          }
        }
        tmpClusterWeights = tmpClusterWeights / length(dataOne)
        eachPointLoglik = rep(0, length(dataOne))
        for (k in 1:length(dataOne)) {
          eachPointLoglik[k] = log(sum(tmpClusterWeights * matrixOfLikeliks[k,,]))
          if (is.na(eachPointLoglik[k])){eachPointLoglik[k] = -10**100}
        }
        
        tmpCurrentLoglik = sum(eachPointLoglik)
        BIC = -2 * tmpCurrentLoglik + log(length(dataOne)) * length(which(tmpClusterWeights != 0))
        if (BIC < oldBIC) {
          currentLoglik = tmpCurrentLoglik
          clusterWeights = tmpClusterWeights
          oldBIC = BIC
        } else {
          break
        }
      }
      
      return(list(currentLoglik, clusterWeights))
    }
  }
}

