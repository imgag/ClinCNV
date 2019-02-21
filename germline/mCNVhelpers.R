
EstimateModeSimple <- function(x) {
  tmpx = x[which(x > 0.5)]
  if (length(tmpx) < 5) {
    return(0)
  }
  
  density_of_x <-  density(tmpx, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  distToMode <- abs(tmpx - mu)
  threshold = max(quantile(distToMode, 0.2), (1 - sqrt(11/12)) / 2)
  if (length(which(distToMode <= threshold) > 40)) {
    return(median(tmpx[which(distToMode <= threshold)]))
  } else {
    lehmanHodges(tmpx[which(distToMode <= threshold)])
  }
}



likelihoodOfGaussianMixture <- function(location, value , sd_to_start, robustPercentage, minSizeOfCluster, esimtatedVarianceFromSampleNoise, lowerBoundSD) {
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
      multiplierForLocation = median(locationTmp / location)
      location = location * multiplierForLocation
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
  if (length(which(newlm1$residuals > 2.5 * QnResid1)) < 0.05 * length(whichBothNonHomo) & length(which(newlm2$residuals > 2.5 * QnResid2)) < 0.05 * length(whichBothNonHomo)) {
    return(T)
  }
  plot(covOne ~ covTwo)
  return(F)
}


findFinalState <- function(coverageNeededToCheck, medianOfCoverage, toyBedFilePolymorphCurrent, multipliersSamples) {
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
  possibleLocations = round(medianOfCoverage ** 2 * 2)
  possibleLocations = unique(c(possibleLocations, possibleLocations + -3:3))
  for (j in 1:length(locations)) {
    if (!j %in% possibleLocations) next
    vecOfMeans = locations[[j]]
    vecOfMeans = vecOfMeans[which(vecOfMeans > min(coverageSummarised[notHomozygousDeletions]) - 0.25 & vecOfMeans < max(coverageSummarised[notHomozygousDeletions]) + 0.25)]
    if (length(vecOfMeans) == 1) next
    likelikAndWeights = likelihoodOfGaussianMixture(vecOfMeans, coverageSummarised[notHomozygousDeletions], sdNormalised, 
                                                    0.05 * (length(notHomozygousDeletions)) / length(coverageSummarised), 0.1, 
                                                    multipliersSamples[notHomozygousDeletions], lowerBoundOfSD)
    firstSignifCluster = which.min(likelikAndWeights[[2]] > 5 / length(notHomozygousDeletions))
    if (length(notHomozygousDeletions) < length(coverageSummarised)) {
      firstSignifCluster = 1
    }
    lastSignifCluster = which.max(likelikAndWeights[[2]] > 5 / length(notHomozygousDeletions))
    if (firstSignifCluster < lastSignifCluster - 1)
      if (min(likelikAndWeights[[2]][(firstSignifCluster + 1):(lastSignifCluster - 1)]) < 1 / length(notHomozygousDeletions)) next
    tmpLikelik =  (
      -2 * likelikAndWeights[[1]] + (length(which(likelikAndWeights[[2]] > 0.1 / length(notHomozygousDeletions)))  + 1) * log(length(notHomozygousDeletions))
    )
    #normalmixEM(coverageSummarised[notHomozygousDeletions], mu=vecOfMeans)
    #print(tmpLikelik)
    #print(j)
    bestWeightCheckForEvenDominance = c(length(coverageSummarised) - length(notHomozygousDeletions), likelikAndWeights[[2]] * length(notHomozygousDeletions)) + 0.1
    bestWeightCheckForEvenDominance = bestWeightCheckForEvenDominance / sum(bestWeightCheckForEvenDominance)
    bestLocEvenDominance = c(0, vecOfMeans)
    bestLocEvenDominance = bestLocEvenDominance ** 2 * j
    if (sum(bestWeightCheckForEvenDominance[which(bestLocEvenDominance %% 2 == 0)]) < 0.4) next

    if (tmpLikelik < bestLikelik | is.null(bestSD)) {
      if (likelikAndWeights[[3]] < 1 - sqrt((j-1)/j)) {
        bestLikelik= tmpLikelik
        bestLoc = vecOfMeans
        bestSD = likelikAndWeights[[3]]
        bestWeight = likelikAndWeights[[2]]
        bestDivisor = j
      }
    }
  }
  bestLoc = unique(c(0, bestLoc))
  bestLoc[which(bestLoc ==0)] = median(coverageSummarised[which(coverageSummarised <= 0.25)])
  bestWeight = c(length(coverageSummarised) - length(notHomozygousDeletions), bestWeight * length(notHomozygousDeletions)) + 1
  bestWeight = bestWeight / sum(bestWeight)
  
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  copy_number_likeliks <- abs(sweep(rep.row(coverageSummarised, length(bestLoc)), 1, bestLoc, FUN="-"))
  copy_number_likeliks <- t(apply(copy_number_likeliks, 1, function(x){dnorm(x, sd = bestSD * multipliersSamples)}))
  copy_number_likeliks <- sweep(copy_number_likeliks, 1, bestWeight, FUN="*")
  copy_number <- bestLoc[apply(copy_number_likeliks, 2, which.max)]
  copy_number = round(copy_number ** 2 * bestDivisor)
  coloursP = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414, 337)]
  coloursP = c(coloursP, coloursP)
  diagnosticPlot = (length(which(copy_number != 2)) > 0.05 * length(copy_number))
  if (diagnosticPlot == T) {
    fileName = paste(toyBedFilePolymorphCurrent[1,1], toyBedFilePolymorphCurrent[1,2], toyBedFilePolymorphCurrent[nrow(toyBedFilePolymorphCurrent),3], sep="_")
    png(paste0(fileName, ".png"), width=length(copy_number) * 3, height=800)
    plot(coverageSummarised, col="black", pch=21,bg=coloursP[(copy_number + 1)])
    abline(h=bestLoc)
    dev.off()
  }
  #plot(density(coverageSummarised, bw='SJ'))
  return(copy_number)
}
