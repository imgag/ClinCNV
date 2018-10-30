fast_dt_list <- function(degreesOfFreedom) {
  values <- seq(from = 0.0, to = 10000.0, by=1)
  vect_of_t_likeliks <- dt(values / 1000, df=degreesOfFreedom)
  return((vect_of_t_likeliks))
}


return_likelik <- function(x) {
  x = as.vector(x)
  x = round(abs(x * 1000)) + 1
  x = replace(x, which(x >= length(vect_of_t_likeliks)), length(vect_of_t_likeliks) - 1)
  return(vect_of_t_likeliks[x])
}


EstimateModeSimple <- function(x) {
  density_of_x <-  density(x, kernel="gaussian")
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
  
  # sample size normalisation
  normalization_factors <- rep(1, ncol(coverages))
  
  
  logarithms <- apply(coverages[autosomes,], 1, function(x) {return(sum(log(x)) / ncol(coverages))})
  root_of_mth_power <- exp(logarithms)
  
  for (i in 1:ncol(coverages)) {
    array_to_calculaste_medians <- coverages[autosomes,i] / root_of_mth_power
    normalization_factors[i] <- median(array_to_calculaste_medians)
  }
  for (i in 1:ncol(coverages)) {
    coverages[,i] <- (coverages[,i] / (normalization_factors [i]))
  }
  
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
        position <- which(startsWith(names(allowedChroms), prefix=tumorName))
        if (length(position) == 1) {
        allowedChromosomesAutosomesOnly = which(!info[,1] %in% c("chrX", "chrY") & info[,1] %in% allowedChroms[[position]])
        } else {
          allowedChromosomesAutosomesOnly = which(!info[,1] %in% c("chrX", "chrY"))
        }
        vector_of_gcs_in_allowed_chroms = intersect(vector_of_gc, allowedChromosomesAutosomesOnly)
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
















find_final_state <- function(start, end, initial_state, matrix_of_likeliks) {
  super_small_likelik = -10^20
  sweeped_matrix <- sweep(matrix_of_likeliks[(start + 1):(end - 1),,drop=F], 1, matrix_of_likeliks[(start + 1):(end - 1),initial_state])
  res_within <- apply(sweeped_matrix, 2, sum)
  cn_state_by_central_points <- which.min(res_within)
  
  only_central=F
  sweeped_matrix <- sweep(matrix_of_likeliks[(start):(end),,drop=F], 1, matrix_of_likeliks[(start):(end),initial_state])
  res <- apply(sweeped_matrix, 2, sum)
  likelihood_score_including_all_tiles <- min(res[cn_state_by_central_points], res_within[cn_state_by_central_points])
  if (likelihood_score_including_all_tiles == res_within[cn_state_by_central_points]) {
    only_central=T
  }
  
  if (!only_central) {
    sweeped_matrix <- sweep(matrix_of_likeliks[(start):(end),,drop=F], 1, matrix_of_likeliks[(start):(end),initial_state])
    res <- apply(sweeped_matrix, 2, sum)
    likelihood_score_all_tiles_read_depth_only <- res[cn_state_by_central_points]
    for_output <- apply((matrix_of_likeliks[(start + 1):(end - 1),,drop=F] / -2) * log10(exp(0)), 2, sum)
  } else {
    sweeped_matrix <- sweep(matrix_of_likeliks[(start + 1):(end - 1),,drop=F], 1, matrix_of_likeliks[(start + 1):(end - 1),initial_state])
    res <- apply(sweeped_matrix, 2, sum)
    likelihood_score_all_tiles_read_depth_only <- res[cn_state_by_central_points]
    for_output <- apply((matrix_of_likeliks[(start + 1):(end - 1),,drop=F] / -2) * log10(exp(0)), 2, sum)
  }

  
  return(c(likelihood_score_including_all_tiles, cn_state_by_central_points, likelihood_score_all_tiles_read_depth_only, for_output))
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

find_one_CNV <- function(j, k, main_state, threshold, matrix_of_likeliks, min_CNV_len) {
  # sweeping the likelihoods
  subset_matrix <- -1 * matrix_of_likeliks[j:k,,drop=F]
  matrix_of_BFs <- sweep(subset_matrix, 1, subset_matrix[,main_state], FUN="-")
  coords_of_CNVs <- c(0,0,0,0)
  best_bf <- 0
  value <- threshold
  sequence_for_iteration = seq(1:ncol(matrix_of_BFs))
  sequence_for_iteration = sequence_for_iteration[-main_state]
  for (i in sequence_for_iteration) {
    res <- maxSubArraySum(matrix_of_BFs[,i])
    res[2] <- j + res[2] - 1
    res[3] <- j + res[3] - 1
    if (res[1] > max(value, best_bf)) {
      best_bf <- res[1]
      coords_of_CNVs <- c(res, i)
    }
  }
  if (coords_of_CNVs != c(0,0,0,0)) {
    return(coords_of_CNVs)
  }
  return(c(0,0,0,0))
}


find_all_CNVs <- function(minimum_length_of_CNV, threshold, price_per_tile, initial_state, matrix_of_likeliks, very_initial_state) {
  vector_of_regions <- matrix(c(-10, 1, nrow(matrix_of_likeliks), initial_state), nrow=1, ncol=4)
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
      found_CNV = find_one_CNV(start, end, current_region_to_look_for_CNVs[4], threshold, matrix_of_likeliks, minimum_length_of_CNV)
      if (found_CNV[4] == 0)
        flag_not_found_or_too_short = T
    } else {
      flag_not_found_or_too_short = T
    }
    if (flag_not_found_or_too_short) {
      # if we do not segment further we add CNV to the list - found CNV is not significant or 
      result_CNV <- current_region_to_look_for_CNVs
      if (current_region_to_look_for_CNVs[4] != initial_state & end - start >= minimum_length_of_CNV){
        bf_and_state <- find_final_state(start, end, very_initial_state, matrix_of_likeliks)
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
      evaluate_segment_further = T
     
        for (k in found_CNV[2]:found_CNV[3]) {
          matrix_of_likeliks[k,] = 0
        }

      matrix_for_calculations <- -matrix_of_likeliks[start:end,] + matrix_of_likeliks[start:end, initial_state]
      if (!is.null(nrow(matrix_for_calculations))) {
        matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
      } 
      determined_state <- which.max(matrix_for_calculations)
      BF <- max(matrix_for_calculations)
      current_region_to_look_for_CNVs <- c(BF, start, end, determined_state)
      if (evaluate_segment_further)
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
        matrix_for_calculations <- -matrix_of_likeliks[start:(start_of_CNV - 1),] + matrix_of_likeliks[start:(start_of_CNV - 1), current_region_to_look_for_CNVs[4]]
        if (!is.null(nrow(matrix_for_calculations))) {
          matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
        } 
        determined_state <- which.max(matrix_for_calculations)
        BF <- max(matrix_for_calculations)
        left_part <- c(BF, start, start_of_CNV - 1, determined_state)
        vector_of_regions <- rbind(vector_of_regions, left_part)
      }
      if (end - end_of_CNV >= minimum_length_of_CNV) { # if right part is big enough to add! CAUTION!!!
        matrix_for_calculations <- -matrix_of_likeliks[start:(start_of_CNV - 1),] + matrix_of_likeliks[start:(start_of_CNV - 1), current_region_to_look_for_CNVs[4]]
        if (!is.null(nrow(matrix_for_calculations))) {
          matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
        } 
        determined_state <- which.max(matrix_for_calculations)
        BF <- max(matrix_for_calculations)
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
  maxCopyNumber = 6
  if (nrow(toyBedFile) == 0 || length(dotsCoords) == 0) {
    return(0)
  } else {
    ### MAKE ANNOTATION TO TRACK FILES
    makeTrackAnnotation <- function(fileName) {
      if (!file.exists(fileName)) {
        file.create(fileName)
        fileConn<-file(fileName)
        writeLines(c("#type=GENE_EXPRESSION",
                     paste0("#track graphtype=points name=\"", ID, "\" color=0,0,255 altColor=255,0,0 maxHeightPixels=80:80:80 viewLimits=0:6 yLineMark=2 yLineOnOff=on"),
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
    copyNumberValues <- round(reverseFunctionUsedToTransform(dotsCoords), digits=2)
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






cleanDatasetFromLowCoveredFiles <- function(normal) {
  medians <- apply(sqrt(normal), 1, median)
  minAllowedCoverage = max(quantile(medians, 0.01), 0.05)
  rowsToRemove <- which(medians < minAllowedCoverage)
  return(rowsToRemove)
}