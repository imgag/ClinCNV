
startsWith <- function(x, prefix) {
  if (substring(x, 1, nchar(prefix)) == prefix) {
    return(T)
  } else {
    return(F)
  }
}

ReadFileFast <- function(fileName, header=T) {
  if (header) {
    colnames <- strsplit(readLines(fileName, n=1), "\t")[[1]]
    setnames(localDf <- as.data.frame(fread(fileName, skip=1, header=F, stringsAsFactors = F, sep="\t")), colnames)
  } else {
    localDf <- as.data.frame(fread(fileName, header=F, stringsAsFactors = F, sep="\t"))
  }
  return(localDf)
}

fast_dt_list <- function(degreesOfFreedom) {
  values <- seq(from = 0.0, to = 10000.0, by=1)
  vect_of_t_likeliks <- dt(values / 1000, df=degreesOfFreedom)
  return((vect_of_t_likeliks))
}


return_likelik <- function(x) {
  x = as.vector(x)
  x = round(abs(x * 1000))
  x[which(x < 1)] = 1
  x = replace(x, which(x >= length(vect_of_norm_likeliks)), length(vect_of_norm_likeliks) - 1)
  return(vect_of_norm_likeliks[x])
}


fast_dnorm_list <- function() {
  values <- seq(from = 0.0, to = 10000.0, by=1)
  vect_of_norm_likeliks <- dnorm(values / 1000)
  return((vect_of_norm_likeliks))
}


return_t_likelik <- function(x) {
  x = as.vector(x)
  x = round(abs(x * 1000)) + 1
  x = replace(x, which(x >= length(vect_of_t_likeliks)), length(vect_of_t_likeliks) - 1)
  return(vect_of_t_likeliks[x])
}

getCytobands <- function(fileName) {
  cytobands <- read.table(fileName,stringsAsFactors = F, header = F, sep="\t", row.names=NULL)
  left_borders <- vector(mode="list", length=nrow(cytobands)/2)
  right_borders <- vector(mode="list", length=nrow(cytobands)/2)
  ends_of_chroms <- vector(mode="list", length=nrow(cytobands)/2)
  odd_numbers <- seq(from=1, to=nrow(cytobands), by=2)
  even_numbers <- seq(from=2, to=nrow(cytobands), by=2)
  names(left_borders) = as.vector(t(cytobands[,1]))[odd_numbers]
  names(right_borders) = as.vector(t(cytobands[,1]))[even_numbers]
  names(ends_of_chroms) = as.vector(t(cytobands[,1]))[even_numbers]
  for (i in 1:length(odd_numbers)) {
    left_borders[[i]] = cytobands[odd_numbers[i], 2]
    right_borders[[i]] = cytobands[even_numbers[i], 3]
    ends_of_chroms[[i]] = cytobands[even_numbers[i], 5]
  }
  return(list(left_borders, right_borders, ends_of_chroms))
}


EstimateModeSimple <- function(x) {
  density_of_x <-  density(x, kernel="gaussian", bw="SJ")
  mu = density_of_x$x[which.max(density_of_x$y)]
  mu
}

lehmanHodges <- function(x) {
  allCombs <- combn(x, 2, FUN=mean)
  averages <- c(x, allCombs)
  return(median(averages))
}


gc_and_sample_size_normalise <- function(info, coverages, averageCoverage=T, allowedChroms=NULL) {
  # allowedChroms is the list, index in the list = column in coverages
  for (j in 1:ncol(coverages)) {
    smallValues <- which(coverages[,j] < 10^-20)
    coverages[smallValues, j] <- (10^-20)
  }
  
  if (!averageCoverage) {
    for (i in 1:nrow(coverages)) {
      coverages[i,] = coverages[i,] / (info[i,3] - info[i,2])
    }
  }
  
  autosomes <- which(!info[,1] %in% c("chrX", "chrY"))
  

  
  coverages = log2(coverages)
  
  gcs <- info[,4]
  uniques_gcs <- unique(gcs)
  uniquesGcsToExclude = c()
  for (i in 1:length(uniques_gcs)) {
    if (length(which(gcs == uniques_gcs[i])) < 100) {
      uniquesGcsToExclude = c(uniquesGcsToExclude, i)
    }
  }
  if (!length(uniquesGcsToExclude) == 0)
  uniques_gcs = uniques_gcs[-uniquesGcsToExclude]
  print(paste("Percentage of regions remained after GC correction:", length(which(gcs %in% uniques_gcs)) / length(gcs)))
  
  if (is.null(allowedChroms)) {
    allowedChromosomesAutosomesOnly = autosomes
    gc_normalisation_factors = foreach (i = 1:length(uniques_gcs), .combine="rbind", .export=c("EstimateModeSimple", "lehmanHodges")) %dopar% {
      gc()
      curr_gc = uniques_gcs[i]
      vector_of_gc <- which(gcs == curr_gc)
      vector_of_gcs_in_allowed_chroms = intersect(vector_of_gc, allowedChromosomesAutosomesOnly)
      if (length(vector_of_gcs_in_allowed_chroms) >= 50) {
        gc_norm_factor <- as.vector(apply(coverages[vector_of_gcs_in_allowed_chroms,], 2, median))
      } else {
        gc_norm_factor = as.vector(apply(coverages[vector_of_gcs_in_allowed_chroms,], 2, lehmanHodges))
      }
      gc_norm_factor
    }
  } else {
    gc_normalisation_factors = foreach (i = 1:length(uniques_gcs), .combine="rbind", .export=c("EstimateModeSimple", "lehmanHodges")) %dopar% {
      
      
      gc()
      curr_gc = uniques_gcs[i]
      vector_of_gc <- which(gcs == curr_gc)
      gc_norm_factor = rep(1, ncol(coverages))
      for (j in 1:ncol(coverages)) {
        tumorName = colnames(coverages)[j]
        position <- which(substring(names(allowedChroms), 1, nchar(tumorName)) == tumorName)
        if (length(position) == 1) {
          allowedChromosomesAutosomesOnly = c()
          for (allowedArm in allowedChroms[[position]]) {
            splittedValue <- strsplit(allowedArm, "-")
            chrom = splittedValue[[1]][1]
            if (!chrom %in% c("chrX", "chrY", "X", "Y")) {
              startOfArm = as.numeric(splittedValue[[1]][2])
              endOfArm = as.numeric(splittedValue[[1]][3])
              allowedChromosomesAutosomesOnly = union(allowedChromosomesAutosomesOnly, which(info[,1] == chrom &
                                                                                               info[,2] >= startOfArm &
                                                                                               info[,3] <= endOfArm))
            }
          }
        } else {
          allowedChromosomesAutosomesOnly = which(!info[,1] %in% c("chrX", "chrY"))
        }
        vector_of_gcs_in_allowed_chroms = intersect(vector_of_gc, allowedChromosomesAutosomesOnly)
        borderOfDistnace = 0.00
        while (length(vector_of_gcs_in_allowed_chroms) < 20) {
          borderOfDistnace = borderOfDistnace + 0.01
          tmpVectorOfGC = which(gcs %in% c(curr_gc + seq(from=-borderOfDistnace, to=borderOfDistnace, by=0.01)))
          vector_of_gcs_in_allowed_chroms = intersect(tmpVectorOfGC, allowedChromosomesAutosomesOnly)
        }
        if (length(vector_of_gcs_in_allowed_chroms) >= 50) {
          gc_norm_factor[j] = median(coverages[vector_of_gcs_in_allowed_chroms,j])
        } else {
          gc_norm_factor[j] = lehmanHodges(coverages[vector_of_gcs_in_allowed_chroms,j])
        }
      }
      gc_norm_factor
    }
  }

  
  for (i in 1:length(uniques_gcs)) {
    factorGC = gc_normalisation_factors[i,]
    rowsWithThatGC <- which(gcs == uniques_gcs[i])
    for (rowCoord in rowsWithThatGC) {
      coverages[rowCoord,] = coverages[rowCoord,] - gc_normalisation_factors[i,]
    }
  }
  
  if (length(which(! gcs %in% uniques_gcs)) > 0) {
    coverages <- coverages[-which(! gcs %in% uniques_gcs),]
    info <- info[-which(! gcs %in% uniques_gcs),]
  }
  
  return(list(2 ** (coverages), info))
}









getBordersOfChromosomes <- function(bedFile) {
  prevChrom = "chrN"
  bordersOfChroms <- c()
  for (j in 1:nrow(bedFile)) {
    if (bedFile[j,1] != prevChrom) {
      bordersOfChroms <- c(bordersOfChroms, j)
      prevChrom = bedFile[j,1]
    }
  }
  return(bordersOfChroms)
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
















find_final_state <- function(start, end, initial_state, matrix_of_likeliks_local, blocked_states, penalties) {
  super_small_likelik = -10^20
  sweeped_matrix <- sweep(matrix_of_likeliks_local[min(start + 1, end):(max(end - 1, start)),,drop=F], 1, matrix_of_likeliks_local[min(start + 1, end):(max(end - 1, start)),initial_state])
  res_within <- apply(sweeped_matrix, 2, sum)
  if (length(penalties) == ncol(sweeped_matrix)) {
    res_within = res_within + penalties
  }
  res_within[blocked_states] = max(res_within) + 1
  cn_state_by_central_points <- which.min(res_within)
  
  only_central=F
  sweeped_matrix <- sweep(matrix_of_likeliks_local[(start):(end),,drop=F], 1, matrix_of_likeliks_local[(start):(end),initial_state])
  res <- apply(sweeped_matrix, 2, sum)
  likelihood_score_including_all_tiles <- min(res[cn_state_by_central_points], res_within[cn_state_by_central_points])
  if (likelihood_score_including_all_tiles == res_within[cn_state_by_central_points]) {
    only_central=T
  }
  
  if (!only_central) {
    sweeped_matrix <- sweep(matrix_of_likeliks_local[(start):(end),,drop=F], 1, matrix_of_likeliks_local[(start):(end),initial_state])
    res <- apply(sweeped_matrix, 2, sum)
    likelihood_score_all_tiles_read_depth_only <- res[cn_state_by_central_points]
    #for_output <- apply((matrix_of_likeliks_local[(start + 1):(end - 1),,drop=F] / -2) * log10(exp(0)), 2, sum)
  } else {
    sweeped_matrix <- sweep(matrix_of_likeliks_local[min(start + 1, end):(max(end - 1, start)),,drop=F], 1, matrix_of_likeliks_local[min(start + 1, end):(max(end - 1, start)),initial_state])
    res <- apply(sweeped_matrix, 2, sum)
    likelihood_score_all_tiles_read_depth_only <- res[cn_state_by_central_points]
    #for_output <- apply((matrix_of_likeliks_local[(start + 1):(end - 1),,drop=F] / -2) * log10(exp(0)), 2, sum)
  }

  
  return(c(likelihood_score_including_all_tiles, cn_state_by_central_points, likelihood_score_all_tiles_read_depth_only))#, for_output))
}

maxSubArraySum <- function(x){
  bestSoFar = 0
  bestNow = 0
  bestStartIndexSoFar = -1
  bestStopIndexSoFar = -1
  bestStartIndexNow = -1
  for (i in 1:length(x)) {
    value = bestNow + x[i]
    if (value > 0) {
      if (bestNow == 0) {
        bestStartIndexNow = i
      }
      bestNow = value
    }
    else
      bestNow = 0
    
    if (bestNow > bestSoFar) {
      bestSoFar = bestNow
      bestStopIndexSoFar = i
      bestStartIndexSoFar = bestStartIndexNow
    }
  }
  return(c(bestSoFar, bestStartIndexSoFar, bestStopIndexSoFar))
}

find_one_CNV <- function(j, k, main_state, threshold, matrix_of_likeliks_local, min_CNV_len, blocked_states) {
  # sweeping the likelihoods
  subset_matrix <- -1 * matrix_of_likeliks_local[j:k,,drop=F]
  matrix_of_BFs <- sweep(subset_matrix, 1, subset_matrix[,main_state], FUN="-")
  coords_of_CNVs <- c(0,0,0,0)
  best_bf <- 0
  value <- threshold
  sequence_for_iteration = seq(1:ncol(matrix_of_BFs))
  sequence_for_iteration = setdiff(sequence_for_iteration, c(blocked_states, main_state))
  detectedCNVs <- foreach (i=sequence_for_iteration, .combine=rbind) %dopar% {
    maxSubArraySum(matrix_of_BFs[,i])
  }
  detectedCNVs = matrix(detectedCNVs, ncol=3)
  resultCNV = detectedCNVs[which.max(detectedCNVs[,1]),]
  if (resultCNV[1] > threshold) {
    resultCNV[2] = j + resultCNV[2] - 1
    resultCNV[3] = j + resultCNV[3] - 1
    resultCNV <- c(resultCNV, sequence_for_iteration[which.max(detectedCNVs[,1])])
  } else {
    resultCNV = coords_of_CNVs
  }
  coords_of_CNVs = unname(resultCNV)
  if (coords_of_CNVs != c(0,0,0,0)) {
    return(coords_of_CNVs)
  }
  return(c(0,0,0,0))
}


find_all_CNVs <- function(minimum_length_of_CNV, threshold, price_per_tile, initial_state, matrix_of_likeliks_local, very_initial_state, blocked_states=c(), penalties=c()) {
  vector_of_regions <- matrix(c(-10, 1, nrow(matrix_of_likeliks_local), initial_state), nrow=1, ncol=4)
  found_CNVs <- matrix(nrow=0, ncol=11)
  i = 1
  counter = 0
  while(i <= nrow(vector_of_regions)){
    current_region_to_look_for_CNVs = vector_of_regions[i,]
    start = current_region_to_look_for_CNVs[2]
    end = current_region_to_look_for_CNVs[3]
    allowed_length = max(3, minimum_length_of_CNV)
    flag_not_found_or_too_short = F
    if (end - start > allowed_length){
      found_CNV = find_one_CNV(start, end, current_region_to_look_for_CNVs[4], threshold, matrix_of_likeliks_local, minimum_length_of_CNV, blocked_states)
      if (found_CNV[4] == 0)
        flag_not_found_or_too_short = T
    } else {
      flag_not_found_or_too_short = T
    }
    if (flag_not_found_or_too_short) {
      # if we do not segment further we add CNV to the list - found CNV is not significant or 
      result_CNV <- current_region_to_look_for_CNVs
      if (current_region_to_look_for_CNVs[4] != initial_state & end - start >= minimum_length_of_CNV){
        bf_and_state <- find_final_state(start, end, very_initial_state, matrix_of_likeliks_local, blocked_states, penalties)
        result_CNV[1] = bf_and_state[1]
        result_CNV[4] = bf_and_state[2]
        likelik_score_read_depth_only <- bf_and_state[3]
        if (bf_and_state[2] != initial_state & bf_and_state[1] < min(-threshold,  -price_per_tile * (end - start + 1) )) {
          if (likelik_score_read_depth_only < -threshold) {
            found_CNVs <- rbind(found_CNVs, result_CNV)
          }
        }
      }
    } else if (found_CNV[4] != 0 & found_CNV[3] - found_CNV[2] < minimum_length_of_CNV) {
      # found CNV is too short!
      factor = found_CNV[1] / (threshold - 1)
      matrix_of_likeliks_local[found_CNV[2]:found_CNV[3],] = as.matrix(matrix_of_likeliks_local[found_CNV[2]:found_CNV[3],] / factor)

      matrix_for_calculations <- sweep(matrix_of_likeliks_local[start:end,,drop=F], 1, matrix_of_likeliks_local[start:end, initial_state, drop=F], FUN="-")
      if (!is.null(nrow(matrix_for_calculations))) {
        matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
      } 
      determined_state <- which.min(matrix_for_calculations)
      BF <- -1 * min(matrix_for_calculations)
      current_region_to_look_for_CNVs <- c(BF, start, end, determined_state)
      vector_of_regions <- rbind(vector_of_regions, current_region_to_look_for_CNVs)
    } else {
      # found CNV is significant and long
      start_of_CNV <- found_CNV[2]
      end_of_CNV <- found_CNV[3]
      new_state <- found_CNV[4]
      # segment found CNV itself
      if (end_of_CNV - start_of_CNV >= minimum_length_of_CNV) {
        vector_of_regions <- rbind(vector_of_regions, found_CNV)
      }
      if (start_of_CNV - start >= minimum_length_of_CNV) { # if left part is big enough to add! CAUTION!!!
        matrix_for_calculations <- sweep(matrix_of_likeliks_local[start:max(start, start_of_CNV - 1),,drop=F], 2, matrix_of_likeliks_local[start:max(start, start_of_CNV - 1), current_region_to_look_for_CNVs[4],drop=F], FUN="-")
        if (!is.null(nrow(matrix_for_calculations))) {
          matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
        } 
        determined_state <- which.min(matrix_for_calculations)
        BF <- -1 * min(matrix_for_calculations)
        left_part <- c(BF, start, start_of_CNV - 1, determined_state)
        vector_of_regions <- rbind(vector_of_regions, left_part)
      }
      if (end - end_of_CNV >= minimum_length_of_CNV) { # if right part is big enough to add! CAUTION!!!
        matrix_for_calculations <- sweep(matrix_of_likeliks_local[min(end, end_of_CNV + 1):end,,drop=F], 2, matrix_of_likeliks_local[min(end, end_of_CNV + 1):end, current_region_to_look_for_CNVs[4],drop=F], FUN="-")
        if (!is.null(nrow(matrix_for_calculations))) {
          matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
        } 
        determined_state <- which.min(matrix_for_calculations)
        BF <- -1 * min(matrix_for_calculations)
        right_part <- c(BF, end_of_CNV + 1, end, determined_state)
        vector_of_regions <- rbind(vector_of_regions, right_part)
      }
    }
    i = i + 1
  }
  return(found_CNVs)
}










outputSegmentsAndDotsFromListOfCNVs <- function(toyBedFile, foundCNVs, startOfChromPiece, endOfChromPiece, outputFileNameCNVs, 
                                         outputFileNameDots, ID, dotsCoords, reverseFunctionUsedToTransform, cn_states) {
  maxCopyNumber = 8
  if (nrow(toyBedFile) == 0 || length(dotsCoords) == 0) {
    return(0)
  } else {
    ### MAKE ANNOTATION TO TRACK FILES
    makeTrackAnnotation <- function(fileName) {
      if (!file.exists(fileName)) {
        file.create(fileName)
        fileConn<-file(fileName)
        writeLines(c("#type=GENE_EXPRESSION",
                     paste0("#track graphtype=points name=\"", ID, "\" color=0,0,255 altColor=255,0,0 maxHeightPixels=80:80:80 viewLimits=0:2:8 yLineMark=2 yLineOnOff=on"),
                     paste("ID", "chr", "start", "end", "CN", "loglik", "value", sep="\t")), fileConn)
        
        close(fileConn)
      }
    }
    
    makeTrackAnnotation(outputFileNameCNVs)
    makeTrackAnnotation(outputFileNameDots)
    
    if (nrow(toyBedFile) != length(dotsCoords)) {
      print("WARNING: Number of rows in bed file is different from number of dots to output into IGV vis file!")
      return(0)
    }
    copyNumberValues <- round(reverseFunctionUsedToTransform(dotsCoords, toyBedFile[,1]), digits=2)
    likelihoods <- rep(0, length(dotsCoords))
    
    chromosome = toyBedFile[1,1]
    if (nrow(foundCNVs) > 0) {
      if (foundCNVs[1,1] != -1000) {
        for (i in 1:nrow(foundCNVs)) {
          elem = foundCNVs[i,]
          dotsWithinCNV <- elem[2]:elem[3]
          #copyNumberValues[dotsWithinCNV] = cn_states[elem[4]]
          likelihoods[dotsWithinCNV] = elem[1]
          valueToDisplay = min(maxCopyNumber, cn_states[elem[4]])
          copyNumberSegment = matrix(c(ID, chromosome, toyBedFile[elem[2], 2], toyBedFile[elem[3], 3], cn_states[elem[4]], elem[1], valueToDisplay), nrow=1, ncol=7)
          write(paste(copyNumberSegment[1,], collapse="\t"), file=outputFileNameCNVs, append=TRUE)
        }
      }
    }
    for (i in 1:length(dotsCoords)) {
      start = toyBedFile[i,2]
      end = toyBedFile[i,3]
      valueToDisplay = min(maxCopyNumber, copyNumberValues[i])
      copyNumberSegment = matrix(c(ID, chromosome, start, end, copyNumberValues[i], likelihoods[i], valueToDisplay), nrow=1, ncol=7)
      write(paste(copyNumberSegment[1,], collapse="\t"), file=outputFileNameDots, append=TRUE)
    }
  }
 
}






cleanDatasetFromLowCoveredFiles <- function(normal, bedFile) {
  medians <- sapply(1:nrow(normal), function(i) {x=normal[i,]; if (bedFile[i,1] %in% c("chrX", "chrY")) return(quantile(x, 0.9)) else {return(quantile(x,0.5))}})
  minAllowedCoverage = 0
  rowsToRemove <- which(medians <= minAllowedCoverage)
  return(rowsToRemove)
}






lengthBasedNormalization = function(coverage, bedFile, allowedChroms="") {
  lengthBed = round(log2(bedFile[,3] - bedFile[,2]), digits=1)
  topPercent = quantile(lengthBed, 0.95)
  lengthBed[which(lengthBed > topPercent)] = topPercent
  bottomPercent = quantile(lengthBed, 0.05)
  lengthBed[which(lengthBed < bottomPercent)] = bottomPercent
  orderOfLengths = order(lengthBed)
  chroms <- bedFile[orderOfLengths, 1]
  lengthBedOrdered = lengthBed[orderOfLengths]
  chromsToRemoveSex = which(!chroms %in% c("chrX", "chrY"))
  
  
  for (j in 1:ncol(coverage)) {
    
    if (allowedChroms == "") {
      chromsToRemove = chromsToRemoveSex
    } else {
      tumorName = colnames(coverage)[j]
      position <- which(startsWith(names(allowedChroms), prefix=tumorName))
      if (length(position) == 1) {
        chromsToRemove = c()
        for (allowedArm in allowedChroms[[position]]) {
          splittedValue <- strsplit(allowedArm, "-")
          chrom = splittedValue[[1]][1]
          if (!chrom %in% c("chrX", "chrY", "X", "Y")) {
            startOfArm = as.numeric(splittedValue[[1]][2])
            endOfArm = as.numeric(splittedValue[[1]][3])
            chromsToRemove = union(chromsToRemove, which(bedFile[,1] == chrom &
                                                           bedFile[,2] >= startOfArm &
                                                           bedFile[,3] <= endOfArm))
          }
        }
      } else {
        chromsToRemove = chromsToRemoveSex
      }
    }
    lengthBedOrderedLocal = lengthBedOrdered[chromsToRemove]
    
    
    
    
    coverageForNormalization = coverage[orderOfLengths,j]
    medians <- rep(1, length(lengthBed))
    if (length(chromsToRemove) > 0.1 * nrow(bedFile)) {
      listOfMedians = list()
      coverageForNormalizationWithoutBadRegions = coverageForNormalization[chromsToRemove]
      zerosInData = which(coverageForNormalizationWithoutBadRegions < 0.001)
      coverageForNormalizationWithoutZeros = coverageForNormalizationWithoutBadRegions
      if (length(zerosInData > 0)) {
        coverageForNormalizationWithoutZeros = coverageForNormalizationWithoutBadRegions[-zerosInData]
        lengthBedOrderedLocal = lengthBedOrderedLocal[-zerosInData]
      }
      lws = lowess(sqrt(coverageForNormalizationWithoutZeros) ~ lengthBedOrderedLocal, f=2/3)
      
      for (i in 1:length(coverageForNormalization)) {
        lengthOfRegions = lengthBed[i]
        if (as.character(lengthOfRegions) %in% names(listOfMedians)) {
          medians[i] = listOfMedians[[as.character(lengthOfRegions)]]
        } else {
          
          closestX = which.min(abs(lws$x - lengthOfRegions))
          listOfMedians[[as.character(lengthOfRegions)]]  = lws$y[closestX]
          medians[i] = listOfMedians[[as.character(lengthOfRegions)]]
        }
      }
      medians[which(medians == 0)] = median(coverageForNormalizationWithoutBadRegions)
    }
    coverage[,j] = coverage[,j] / (medians ** 2)
  }
  return(coverage)
}



cutX <- function(vec) {
  newVec <- vec
  for (i in 1:length(vec)) {
    elem = vec[i]
    if (startsWith(elem, "X")) {
      newVec[i] = substr(elem, 2, nchar(elem))
    }
  }
  return(newVec)
}



robust_correlation <- function(robust_std, estimation_of_center_x, estimation_of_center_y, x, y) {
  square_root_of_two <- sqrt(2) 
  std_of_x <- robust_std(x)
  std_of_y <- robust_std(y)
  if (std_of_x == 0 | std_of_y == 0) {return(0)}
  first_component = (x - estimation_of_center_x) / (square_root_of_two * std_of_x)
  second_component = (y - estimation_of_center_y) / (square_root_of_two * std_of_y)
  u = first_component + second_component
  v = first_component - second_component
  var_of_u = robust_std(u) ** 2
  var_of_v = robust_std(v) ** 2
  r = (var_of_u - var_of_v) / (var_of_u + var_of_v + 10**-10)
  return(r)
}

robust_correlation_short <- function(robust_std_x, robust_std_y, x, y, robust_std) {
  # use only if data is centered around 0!
  square_root_of_two <- sqrt(2) 
  std_of_x <- robust_std_x
  std_of_y <- robust_std_y
  if (std_of_x == 0 | std_of_y == 0) {return(0)}
  first_component = (x) / (square_root_of_two * std_of_x)
  second_component = (y) / (square_root_of_two * std_of_y)
  u = first_component + second_component
  v = first_component - second_component
  var_of_u = robust_std(u) ** 2
  var_of_v = robust_std(v) ** 2
  r = (var_of_u - var_of_v) / (var_of_u + var_of_v + 10**-10)
  return(r)
}

correlationMatrixForPairedLikelik <- function(x, y, robust_std_x=NULL, robust_std_y=NULL) {
  if (is.null(robust_std_x) & is.null(robust_std_y)) {
    correlationBetweenValues <- robust_correlation(Qn, 0, 0, x, y)
  } else {
    correlationBetweenValues=robust_correlation_short(robust_std_x, robust_std_y, x, y, Qn)
  }
  return(matrix(c(1,correlationBetweenValues,correlationBetweenValues,1), nrow=2))
}


findRegionsToFilerOutDueSystematicallyLowCoverage <- function(normal, tumor=NULL) {
  quantilesOfRowsNorm <- apply(normal, 1, function(x) {quantile(x, 0.9)})
  if (!is.null(tumor)) {
    quantilesOfRowsTum <- apply(tumor, 1, function(x) {quantile(x, 0.9)})
    minimums <- apply(rbind(quantilesOfRowsNorm, quantilesOfRowsTum), 2, min)
    return(which(minimums < 0.3))
  } else {
    return(which(quantilesOfRowsNorm < 0.3))
  }
}


trimValues <- function(values, perc) {
  if (perc > 0.5) perc = 1 - perc
  values[which(values > quantile(values, 1 - perc))] = quantile(values, 1 - perc)
  values[which(values < quantile(values, perc))] = quantile(values, perc)
  return(values)
}


checkForDuplicatesAndRemove <- function(matr, sampleNameToRemain=NULL) {
  summaries <- apply(matr, 2, summary)
  indicesToRemove = c()
  for (i in 1:ncol(summaries)) {
    for (j in 1:ncol(summaries)) {
      if (i != j & identical(summaries[,i], summaries[,j])) {
        if (!is.null(sampleNameToRemain)) {
          if (colnames(matr)[i] == sampleNameToRemain & !colnames(matr)[j] == sampleNameToRemain) {
            indicesToRemove = c(indicesToRemove, j)
          }
          if (colnames(matr)[j] == sampleNameToRemain & !colnames(matr)[i] == sampleNameToRemain) {
            indicesToRemove = c(indicesToRemove, i)
          }
          if  (!colnames(matr)[j] == sampleNameToRemain & !colnames(matr)[i] == sampleNameToRemain) {
            if (j > i) {
              indicesToRemove = c(indicesToRemove, j)
            }
          }
        } else {
          if (j > i) {
            indicesToRemove = c(indicesToRemove, j)
          }
        }
      }
    }
  }
  if (length(indicesToRemove) > 0) {
    print("Attention - duplicates in your dataset!")
    print(paste("We removed duplicates with names", paste(colnames(matr)[unique(indicesToRemove)], collapse = ", ")))
    return(matr[,-indicesToRemove])
  }
 return(matr)
}



determineAverageDepth <- function(rawCoverage, bedFile) {
  avgDepth <- apply(as.matrix(rawCoverage[which(!bedFile[,1] %in% c("chrX","chrY")),]), 2, median)
  avgDepth
}


writeOutLevelOfNoiseVersusCoverage <- function(avgDepth, gcNormalisedCov, bedFile, nameForOutputFile) {
  namesOfOutputFile = colnames(gcNormalisedCov)
  noises <- round(apply(gcNormalisedCov[which(!bedFile[,1] %in% c("chrX","chrY")),], 2, mad), digits = 2)
  tableForOutput <- cbind(namesOfOutputFile, noises, avgDepth)
  write.table(file = nameForOutputFile, tableForOutput, row.names = F, quote = F, sep="\t")
}
