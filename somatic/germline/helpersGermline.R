EstimateModeSimple <- function(x, chrom="", genders=NULL) {
  if (chrom != "") {
    if (chrom %in% c("X", "Y", "chrX", "chrY")) {
      if (is.null(genders)) {
      x = x[which(x > median(x))]
      } else {
        if (chrom == "chrX") {
          x = x[which(genders == "F")]
          if (length(x) < 2) {
            x = x[which(x > median(x))]
          }
        }
        if (chrom == "chrY") {
          x = sqrt(2) * x[which(genders == "M")]
          if (length(x) < 2) {
            return(1)
          }
        }
      }
    }
  }
  if (length(x) > 50) {
    density_of_x <-  density(x, kernel="gaussian")
    mu = density_of_x$x[which.max(density_of_x$y)]
  } else if (length(x) > 30 ){
    mu = median(x)
  } else {
    mu = lehmanHodges(x)
  }
  if (mu < 0.3) {
    mu = median(x)
  }
  mu
}

    


determineSDsOfGermlineSample <- function(x) {
  return(Qn(x))
}

determineSDsOfGermlineProbe <- function(x, i) {
  if (bedFile[i,1] %in% c("chrX", "chrY")) {
    x = x[which(x > median(x))]
  }
  return(Qn(x))
}

esimtateVarianceFromSampleNoise <- function(vecOfSDsParam, numberOfRepetitions) {
  sds <- rep(0, numberOfRepetitions)
  multiplicatorForSN <- rep(0, numberOfRepetitions)
  for (j in 1:numberOfRepetitions) {
    generatedRandomSample <- rnorm(length(vecOfSDsParam), mean=0, sd=vecOfSDsParam)
    sds[j] = sd(generatedRandomSample)
    multiplicatorForSN[j] = sds[j] / Qn(generatedRandomSample)
  }
  finalSDs <- median(sds)
  return(list(median(multiplicatorForSN), vecOfSDsParam / finalSDs))
}



form_matrix_of_likeliks_one_sample <- function(i, j, k, sds, resid, cn_states) {
  vector_of_values <- resid[,k]
  
  vector_of_states <- cn_states
  matrix_of_BFs <- matrix(0, nrow=(j - i + 1), ncol=length(vector_of_states))
  start <- 1
  end <- j - i + 1
  
  matrix_of_BFs = sapply(1:ncol(matrix_of_BFs), function(l) {
    value = return_likelik((vector_of_values - vector_of_states[l]) / (sds) ) / (sds) + 10^-100
    return(-2 * log(value))
  })
  return(matrix_of_BFs)
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



plotFoundCNVs <- function(found_CNVs, toyLogFoldChange, toyBedFile, outputFolder, chrom, cn_states, toySizesOfPointsFromLocalSds, alleleFrequency, plottingOfPNGs) {
  vector_of_states <- cn_states
  cnvsToOutput <- matrix(0, nrow=0, ncol=9)
  if (nrow(found_CNVs) > 0) {
    for (s in 1:nrow(found_CNVs)) {
      CNV_name <- paste(chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], "CN:", vector_of_states[found_CNVs[s,4]], "-2ln(loglik):", found_CNVs[s,1], 
                        ", number of regions:", found_CNVs[s,3] - found_CNVs[s,2] + 1 ,
                        ", length:", toyBedFile[found_CNVs[s,3],3] - toyBedFile[found_CNVs[s,2],2])
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
                           vector_of_states[found_CNVs[s,4]], -1 * round(found_CNVs[s,1]), 
                           found_CNVs[s,3] - found_CNVs[s,2] + 1,
                           format((toyBedFile[found_CNVs[s,3],3] - toyBedFile[found_CNVs[s,2],2]) / 1000, nsmall=3),
                           format(round(alleleFrequency[s], digits=3), nsmall=3),
                           annotationGenes), nrow=1)
      cnvsToOutput = rbind(cnvsToOutput, CNVtoOut)
      
      length_of_repr <- 1000
      

      st <- found_CNVs[s,2]
      fn <- found_CNVs[s,3]
      
      pr = plottingOfPNGs
      if (pr) {
        plot_st <- max(1,st - length_of_repr)
        plot_fn <- min(length(toyLogFoldChange), fn + length_of_repr)
        png(filename=paste0(outputFolder, "/", paste0(CNV_name_to_write, ".png")), type = "cairo", width = 1536, height = 640)
		#bitmap(filename=paste0(outputFolder, "/", paste0(CNV_name_to_write, ".png")) ,width = 1024, height = 640, units = "px" )

        

        plot(toyLogFoldChange[plot_st:plot_fn], main=CNV_name, ylab="Copy Number", xlab=(paste("CNV within Chromosome Arm" )),
             ylim=c(0, 3), cex=toySizesOfPointsFromLocalSds[plot_st:plot_fn], yaxt='n')
        
        axis(side = 2, at = sqrt(cn_states/2), labels = cn_states)
        abline(v=c(st - plot_st, st - plot_st + fn - st + 1), col="red")
        #text(x = (st - 2 * plot_st + fn) / 2, y = (max(toyLogFoldChange[st:fn]) + 0.1), pos=ifelse(found_CNVs[s,4] > 2, 3, 1), 
        #     labels=c(paste("Genes affected:", paste(unique(bedFile[st:fn, 5]), collapse=", " ))), cex=2,
        #     family = "mono")
        #legend(x = (st - 2 * plot_st + fn) / 2, y = (max(toyLogFoldChange[st:fn]) + 0.1), 
        #       legend=c(paste("Genes affected:", paste(unique(bedFile[st:fn, 5]), collapse=", " ))), cex=2)
        
        
        abline(h=sqrt(cn_states/2),lty=2,col=colours,lwd=3)

        points((st - plot_st):(st - plot_st + fn - st) + 1, toyLogFoldChange[st:fn],col="black", pch=21,bg=colours[found_CNVs[s,4]], cex=toySizesOfPointsFromLocalSds[found_CNVs[s,2]:found_CNVs[s,3]])
        
        
        ### EACH POINTS WITH DISTANCE > 10 MB WILL BE SEPARATED BY VERTICAL LINE
        distanceBetweenPoints = 10 ** 6
        for (i in 2:nrow(toyBedFile)) {
          if (toyBedFile[i,2] - toyBedFile[i - 1,2] > distanceBetweenPoints) {
            abline(v = i - 0.5 - plot_st, lty=2, col="grey")
          }
        }
        
        minimumDist <- min(sqrt(cn_states/2) - max(toyLogFoldChange[max(st - 50, plot_st):min(fn + 50, plot_fn)]))
        closestLine <- which(sqrt(cn_states/2) - max(toyLogFoldChange[max(st - 50, plot_st):min(fn + 50, plot_fn)]) == minimumDist)
        text(x = (st - 2 * plot_st + fn) / 2, y = (sqrt(cn_states/2)[closestLine] + 0.1), pos=3, 
             labels=c(paste("Genes affected:", paste(unique(bedFile[st:fn, 5]), collapse=", " ))), cex=1.5,
             family = "mono")
        dev.off()
      }

    }
  }
  return(cnvsToOutput)
}




Determine.gender <- function(normalized.coverage.corrected.gc, probes, clusterNo) {
  set.seed(100)
  if (length(which(probes[,1] == "chrX")) > 100 & length(which(probes[,1] == "chrY")) > 10) {
    matrix_of_values <- cbind(apply(normalized.coverage.corrected.gc[which(probes[,1] == "chrY"),], 2, EstimateModeSimple), apply(normalized.coverage.corrected.gc[which(probes[,1] == "chrX"),], 2, EstimateModeSimple))
    clKmeans <- NULL
    clKmeans <- tryCatch({kmeans(matrix_of_values, centers=matrix(c(0,1,sqrt(1/2), sqrt(1/2)), nrow=2))}, 
                   error = function(e) {
                     NULL
                   })
    if (!is.null(clKmeans)) {
      clusters <- clKmeans$cluster
      clusters[clusters == 1] <- "F"
      clusters[clusters == 2] <- "M"
      png(filename=paste0(opt$out, paste0("/", clusterNo, "genders.png")), width=800, height=800)
      plot(matrix_of_values, col = clKmeans$cluster, xlab="Y Chromsome", ylab="X Chromosome", pch=19, cex=2)
      points(clKmeans$centers, col = 1:2, pch = 8, cex = 10)
      dev.off()
    } else {
      clusters <- rep(2, nrow(matrix_of_values))
      clusters[which(matrix_of_values[,2] < 0.25)] = 1
      clusters[clusters == 1] <- "F"
      clusters[clusters == 2] <- "M"
    }
  } else {
    clusters = rep("M", ncol(normalized.coverage.corrected.gc))
  }
  
  return(clusters)
}







FindRobustMeanAndStandardDeviation <- function(x, genders, chrom, modeEstimated = NA) {
  genders = genders[which(x > 0.1)]
  x = x[which(x > 0.1)]
  if (length(x) < 50) {
    if (chrom != "") {
      if (chrom %in% c("X", "Y", "chrX", "chrY")) {
        if (is.null(genders)) {
          x = x[which(x > median(x))]
        } else {
          if (chrom == "chrX") {
            if (length(which(genders == "F")) <= 5) {
              x = x[which(genders == "M")]
              x = sqrt(2) * x
            } else {
              x = x[which(genders == "F")]
            }
            if (length(x) <= 5) {
              return(matrix(c(1, 0.5, "sd"), nrow=1))
            }
          }
          if (chrom == "chrY") {
            if (length(which(genders == "M")) <= 5) {
              return(matrix(c(sqrt(1/2), 0.5, "sd"), nrow=1))
            } else {
              x = x[which(genders == "M")]
            }
          }
        }
      }
    }

    if (length(x) > 30 ){
      mu = median(x)
    } else {
      mu = lehmanHodges(x)
    }
    if (mu < 0.3) {
      mu = median(x)
    }
    return(matrix(c(mu, Qn(x), "Qn"), nrow=1))
  }

  if (chrom == "chrX") {
    if (length(which(genders == "F")) <= 5) {
      x = x[which(genders == "M")]
      x = sqrt(2) * x
    } else {
      x = x[which(genders == "F")]
    }
  }
  if (chrom == "chrY") {
    
    if (length(which(genders == "M")) <= 5) {
      return(matrix(c(sqrt(1/2), 0.5, "sd"), nrow=1))
    } else {
      x = x[which(genders == "M")]
    }
  }
  

  forMeanX = x
  return(matrix(c(median(forMeanX), Qn(forMeanX), "Qn"), nrow=1))
  
  ### THEN IT MAY BE CHANGED
  if (length(forMeanX) < 10) {
    return(matrix(c(median(forMeanX), Qn(forMeanX), "Qn"), nrow=1))
  }
  if (length(x) < 100) {
    forMeanX = apply(combn(x, 2), 2, mean)
    density_of_x <-  density(forMeanX, bw="SJ", kernel="gaussian")
  } else {
    density_of_x = density(forMeanX, bw=bandwidth, kernel="gaussian")
  }
  if (is.na(modeEstimated)) {
    mu = density_of_x$x[which.max(density_of_x$y)]
  } else {
    mu = modeEstimated
  }
    if (length(x) < 100) {
      density_of_x <-  density(x, bw=bandwidth, kernel="gaussian")
    }

  closest_to_mu <- which.min(abs(density_of_x$x - mu))
  which_are_bigger <- which(density_of_x$y > density_of_x$y[closest_to_mu])
  density_of_x <- as.data.frame(cbind(density_of_x$x, density_of_x$y))
  colnames(density_of_x) <- c("x","y")
  density_of_x[which_are_bigger,] <- density_of_x$y[closest_to_mu]
  EF = max(density_of_x$y)
  
  lower_bound = min(density_of_x$x)
  lower_bound_differs = F
  bounded_on_lower_copy_nuber = which(density_of_x$x < sqrt(mu ** 2 - 1/8))
  if (length(bounded_on_lower_copy_nuber) > 0) {
    start_to_the_left <- max(bounded_on_lower_copy_nuber)
    
    previous_value = density_of_x$y[start_to_the_left]
    for (i in seq(from = start_to_the_left, to=1, by=-1)) {
      AB =  density_of_x$y[i]
      if ((AB > previous_value + 10**-10) | AB < 10**-10) {
        lower_bound = density_of_x$x[i]
        lower_bound_differs = T
        break
      } else {
        previous_value = AB
      }
    }
  } else {
    lower_bound = sqrt(mu ** 2 + 1/8)
    lower_bound_differs = T
  }
  
  upper_bound = max(density_of_x$x)
  upper_bound_differs = F
  bounded_on_higher_copy_nuber = which(density_of_x$x > sqrt(mu ** 2 + 1/16))
  if (length(bounded_on_higher_copy_nuber) > 0) {
    start_to_the_right <- min(bounded_on_higher_copy_nuber)
    previous_value = density_of_x$y[start_to_the_right]
    for (i in seq(from = start_to_the_right, to=length(density_of_x$x), by=1)) {
      AB =  density_of_x$y[i]
      if ((AB > previous_value + 10**-10 | AB < 10**-10)) {
        upper_bound = density_of_x$x[i]
        upper_bound_differs = T
        break
      } else {
        previous_value = AB
      }
    }
  } else {
    upper_bound = sqrt(mu ** 2 + 1/16)
    upper_bound_differs = T
  }
  
  if (upper_bound_differs & lower_bound_differs) {
    dtnorm0 <- function(X, mean, sd, log = TRUE) {msm::dtnorm(X, mean, sd, lower=lower_bound, upper=upper_bound,
                                                         log=T)}
  } else if (!upper_bound_differs & !lower_bound_differs) {
    return(matrix(c(median(x), Qn(x), "Qn"), nrow=1))
  } else if (upper_bound_differs) {
    dtnorm0 <- function(X, mean, sd, log = FALSE) {msm::dtnorm(X, mean, sd, lower = -10**10, upper=upper_bound,
                                                          log=T)}
  } else {
    dtnorm0 <- function(X, mean, sd, log = FALSE) {msm::dtnorm(X, mean, sd, lower=lower_bound, upper=10**10,
                                                          log=T)}
  }
  QnX <- Qn(x)
  data = x[which(x >= lower_bound & x <= upper_bound)]
  if (is.na(modeEstimated)) {
    
    result <- tryCatch({nres=fitdistr(data, dtnorm0, start=list(mean=mean(data), sd=sd(data))); (nres$estimate)}
                       , error = function(e) {return(matrix(c(mu, QnX, "Qn"), nrow=1))})
  } else {
    result <- tryCatch({nres=fitdistr(data, dtnorm0, fix.arg=list(mean=modeEstimated), start=list(mean=modeEstimated, sd=sd(data))); (nres$estimate)}
                       , error = function(e) {return(matrix(c(mu, QnX, "Qn"), nrow=1))})
  }
  if ((as.numeric(result[1]) / mu) > 1.1 | (as.numeric(result[1]) / as.numeric(mu)) < 0.9 | is.na(result[2])) {
    return(matrix(c(mu, QnX, "Qn"), nrow=1))
  }
  # Sometimes we miss one cluster and that's cause to increase of standard deviation
  result[2] = min(QnX, result[2], sd(x))
  if (result[2] == QnX) {
    result = matrix(c(result[1], result[2], "Qn"), nrow=1)
  } else {
    result = matrix(c(result[1], result[2], "sd"), nrow=1)
  }
  result
}








returnClustering <- function(minNumOfElemsInCluster) {
  set.seed(100)
  clustering = rep(0, ncol(normal))
  if (ncol(normal) < 3 * minNumOfElemsInCluster) {
    return(clustering)
  }
  numOfElementsInCluster = minNumOfElemsInCluster
  

  coverageForClustering = sqrt(normal[which(!bedFile[,1] %in% c("chrX","chrY")),])
  sdsOfRegions <- apply(coverageForClustering, 1, sd)
  potentiallyPolymorphicRegions <- which(sdsOfRegions > quantile(sdsOfRegions, 0.75) | sdsOfRegions == 0)
  
  coverageForClustering = (coverageForClustering[-potentiallyPolymorphicRegions,])
  
  if (!is.null(opt$triosFile)) {
    samplesActuallyPlayingRole = c()
    for (trioRow in 1:nrow(trios)) {
      child_number <- which(colnames(coverageForClustering) == trios[trioRow, 1])
      mother_number  <- which(colnames(coverageForClustering) == trios[trioRow, 2])
      father_number  <- which(colnames(coverageForClustering) == trios[trioRow, 3])
      
      sample_name = paste(trios[trioRow,], collapse="-")
      
      if (length(child_number) == 0 | length(father_number) == 0 | length(mother_number) == 0) {
        next()
      } else {
        samplesActuallyPlayingRole = c(samplesActuallyPlayingRole, c(child_number, mother_number, father_number))
        coverageToReplace = apply(coverageForClustering[,c(child_number, mother_number, father_number)], 1, median)
        coverageForClustering[,child_number] = coverageToReplace
        coverageForClustering[,mother_number] = coverageToReplace
        coverageForClustering[,father_number] = coverageToReplace
      }
    }
    n = 3
  } else {
    n = 3
  }

  coverageForClustering = (apply(coverageForClustering, 2, function(x) tapply(x, ceiling(seq_along(x) / n), median)))
  
  corMatrix <- cor(coverageForClustering)
  
  distMatrix <- sqrt(1 - corMatrix)
  hc <- hclust(as.dist(distMatrix), method="ward.D")
  
  numOfClusters = 1
  for (numOfClusters in 2:100) {
    memb <- cutree(hc, k=numOfClusters)
    numOfObservationsInClusters <- table(memb)
    if (sum(numOfObservationsInClusters[numOfObservationsInClusters >= numOfElementsInCluster]) < 0.8 * sum(numOfObservationsInClusters)) {
      break
    }
  }
  
  plot(hc)
  rect.hclust(hc, k=numOfClusters - 1, border="red")
  memb <- cutree(hc, k=numOfClusters - 1)
  numOfObservationsInClusters <- table(memb)
  clustering = memb
  
  significantClusters = which(numOfObservationsInClusters >= numOfElementsInCluster)
  
  for (i in 1:length(numOfObservationsInClusters)) {
    if (numOfObservationsInClusters[i] < numOfElementsInCluster) {
      for (elem in which(clustering == i)) {
        minDist = 10**1000
        closestCluster = significantClusters[1]
        for (signCluster in significantClusters) {
          distanceToSignCluster = mean(distMatrix[which(memb == signCluster),elem])
          if (distanceToSignCluster < minDist) {
            minDist = distanceToSignCluster
            closestCluster = signCluster
          }
        }
        clustering[elem] = closestCluster
      }
    }
  }
  if (!is.null(opt$triosFile)) {
    clustering[-samplesActuallyPlayingRole] = -1
  }
  
  fit <- cmdscale(dist(t(sqrt(normal[-union(potentiallyPolymorphicRegions, which(bedFile[,1] %in% c("chrX","chrY"))),]))),eig=TRUE, k=2) # k is the number of dim

  # plot solution 
  x <- fit$points[,1]
  y <- fit$points[,2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
       main="Metric MDS", type="n")
  text(x, y, labels = row.names(distMatrix), cex=.7, col=clustering)
  return(clustering)
}


