EstimateModeSimple <- function(x, chrom) {
  if (chrom %in% c("X", "Y", "chrX", "chrY")) {
    x = x[which(x > median(x))]
  }
  density_of_x <-  density(x, kernel="gaussian")
  mu = density_of_x$x[which.max(density_of_x$y)]
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



plotFoundCNVs <- function(found_CNVs, toyLogFoldChange, toyBedFile, outputFolder, chrom, cn_states, toySizesOfPointsFromLocalSds) {
  vector_of_states <- cn_states
  cnvsToOutput <- matrix(0, nrow=0, ncol=6)
  if (nrow(found_CNVs) > 0) {
    for (s in 1:nrow(found_CNVs)) {
      CNV_name <- paste(chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], "CN:", vector_of_states[found_CNVs[s,4]], "-2ln(loglik):", found_CNVs[s,1])
      CNV_name_to_write <- paste(colnames(toyLogFoldChange)[sam_no],  chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], "CN",vector_of_states[found_CNVs[s,4]], sep="_")

      
      vectorOfGeneNames = c()
      genesThatHasToBeSeparated = unique(toyBedFile[found_CNVs[s,2]:found_CNVs[s,3],5])
      print(genesThatHasToBeSeparated)
      for (i in 1:length(genesThatHasToBeSeparated)) {
        vectorOfGeneNames = c(vectorOfGeneNames, unlist(strsplit(genesThatHasToBeSeparated[i], split=",")))
      }
      vectorOfGeneNamesTrimmed = c()
      for (elem in vectorOfGeneNames) {
        vectorOfGeneNamesTrimmed = c(vectorOfGeneNamesTrimmed,trimws(elem) )
      }
      annotationGenes <- paste(unique(vectorOfGeneNamesTrimmed), collapse=",")
      print(annotationGenes)
      CNVtoOut <- matrix(c(chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], vector_of_states[found_CNVs[s,4]], -1 * found_CNVs[s,1], annotationGenes), nrow=1)
      print(CNVtoOut)
      cnvsToOutput = rbind(cnvsToOutput, CNVtoOut)
      
      length_of_repr <- 500
      

      st <- found_CNVs[s,2]
      fn <- found_CNVs[s,3]
      
      pr = T
      if (pr) {
        png(filename=paste0(outputFolder, "/", paste0(CNV_name_to_write, ".png")), type = "cairo", width = 1024, height = 640)
		#bitmap(filename=paste0(outputFolder, "/", paste0(CNV_name_to_write, ".png")) ,width = 1024, height = 640, units = "px" )
        print(CNV_name_to_write)
        print(paste0(outputFolder, CNV_name_to_write))
        

        plot(toyLogFoldChange, main=CNV_name, ylab="Copy Number", xlab=(paste("CNV within Chromosome Arm" )),
             ylim=c(0, 3), cex=toySizesOfPointsFromLocalSds, yaxt='n')
        
        axis(side = 2, at = sqrt(cn_states/2), labels = cn_states)
        abline(v=c(found_CNVs[s,2], found_CNVs[s,3]), col="red")
        
        
        
        abline(h=sqrt(cn_states/2),lty=2,col=colours,lwd=3)
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