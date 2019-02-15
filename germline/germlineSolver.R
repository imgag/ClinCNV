library(party)

no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
cl<-makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

cn_states <- 0:8

#load("/Users/gdemidov/Downloads/prepared.RData")
#opt$folderWithScript = "/Users/gdemidov/Tuebingen/clinCNV_dev_new/ClinCNV/"
#opt$out = "/Users/gdemidov/Tuebingen/CLL/tmpResults"
#locationsShiftedLogFoldChanges <- sweep(matrixOfLogFold, 1, locations)


vect_of_t_likeliks <- fast_dt_list(ncol(coverage.normalised) - 1)
vect_of_norm_likeliks <- fast_dnorm_list()
setwd(opt$folderWithScript)





startCoordOfNonInterruptedSegment = 1
shiftsOfCoverage <- c()
colours = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414,30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414)]

folder_name <- paste0(opt$out, "/normal/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}
covar = T

print(paste("We start to estimate covariances between neighboring regions in germline data - may take some time", Sys.time()))
setwd(opt$out)
covarianceTree = returnTreeForCorrelation(coverage.normalised, 
                                          sdsOfGermlineSamples, sdsOfProbes, 
                                          bedFileFiltered)
if (is.null(covarianceTree)) {
  covar = F
}
setwd(opt$folderWithScript)
print(paste("Tree of covariances (using 2 predictors - sum of regions' lengths and log2 of distance between regions) plotted in", opt$out, Sys.time()))


print(paste("Calling started", Sys.time()))
for (sam_no in 1:ncol(coverage.normalised)) {
  sample_name <- colnames(coverage.normalised)[sam_no]
  if (frameworkOff == "offtargetGermline") { 
    if (sample_name %in% colnames(coverage.normalised.off)) {
    sam_no_off = which(colnames(coverage.normalised.off) == sample_name)
    } else {
      sam_no_off = F
    }
  } else {
    sam_no_off = F
  }
  
  if (!is.null(opt$normalSample)) {
    if (!sample_name == opt$normalSample) {
      next
    }
  }
  print(paste("Working with germline sample", sample_name, Sys.time()))
  
  threshold = opt$scoreG
  minimum_length_of_CNV = opt$lengthG
  price_per_tile = 1
  main_initial_state <- 3
  
  
  localSds = sdsOfProbes * sdsOfGermlineSamples[sam_no]
  localSds[which(localSds == 0)] = median(localSds)
  
  if (sam_no_off) {
    localSdsOff = sdsOfProbesOff * sdsOfGermlineSamplesOff[sam_no_off]
    localSdsOff[which(localSdsOff == 0)] = median(localSdsOff)
  }
  
  
  if (!dir.exists(paste0(folder_name, sample_name))) {
    dir.create(paste0(folder_name, sample_name))
  }
  setwd(paste0(folder_name, sample_name))
  
  
  dict_to_output = c()
  
  
  matrix_of_likeliks_for_FDR <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised), sam_no, localSds, coverage.normalised, sqrt(cn_states / 2))
  matrix_of_likeliks <- matrix_of_likeliks_for_FDR
  
  #if (covar) {
  #  listForLikeliks = form_matrix_of_likeliks_one_sample_with_cov(1, ncol(coverage.normalised), sam_no, localSds, coverage.normalised, sqrt(cn_states / 2), covarianceTree, bedFileFiltered, threshold)
  #   matrix_of_likeliks_with_covar <- listForLikeliks[[1]]
  #   bedFileFilteredWithArtificialProbes = listForLikeliks[[2]]
  #   matrix_of_likeliks <- matrix_of_likeliks_for_FDR
  #} else {
  #  matrix_of_likeliks <- matrix_of_likeliks_for_FDR
  #}
  
  if (sam_no_off) {
    matrix_of_likeliks_off = form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised.off), sam_no_off, localSdsOff, coverage.normalised.off, sqrt(cn_states / 2))
    globalBed = rbind(bedFileFiltered, bedFileFilteredOfftarget)
    orderOfBed = order(globalBed[,1], as.numeric(globalBed[,2]))
    globalBed = globalBed[orderOfBed,]
    globalMatrixOfLikeliks = rbind(matrix_of_likeliks, matrix_of_likeliks_off)
    globalMatrixOfLikeliks = globalMatrixOfLikeliks[orderOfBed,]
  } else {
    globalBed = bedFileFiltered
    globalMatrixOfLikeliks = matrix_of_likeliks
  }

  sizesOfPointsFromLocalSds <- 0.1 / localSds 
  
  
  
  numberOfCNVsIsSufficientlySmall = F
  iterations = 0
  maxIteration = opt$maxNumIter
  vectorWithNumberOfOutliers <- c()
  vectorOfZScores = rep(0, nrow(bedFileFiltered))
  vectorOfZScores[which(!bedFileFiltered[,1] %in% c("chrX","chrY"))] <- (coverage.normalised[which(!bedFileFiltered[,1] %in% c("chrX","chrY")),sam_no] - 1) / localSds
  if (genderOfSamples[sam_no] == "F") {
    vectorOfZScores[which(bedFileFiltered[,1] %in% c("chrX"))] <- c(vectorOfZScores, (coverage.normalised[which(bedFileFiltered[,1] %in% c("chrX")),sam_no] - 1) / localSds)
  } else {
    vectorOfZScores[which(bedFileFiltered[,1] %in% c("chrX","chrY"))] <- c(vectorOfZScores, (coverage.normalised[which(bedFileFiltered[,1] %in% c("chrX","chrY")),sam_no] - sqrt(1/2)) / localSds)
  }
  while (!numberOfCNVsIsSufficientlySmall & iterations < maxIteration) {
    initial_state = main_initial_state
    found_CNVs_total <- matrix(0, nrow=0, ncol=9)
    colnames(found_CNVs_total) <- c("#chr", "start", "end", "CN_change", "loglikelihood", "no_of_regions", "length_KB", "potential_AF", "genes")

    iterations = iterations + 1
    for (l in 1:length(left_borders)) {
      chrom = names(left_borders)[l]
      if (chrom == "chrX" & genderOfSamples[sam_no] == "M") {
        initial_state <- 2
      } else if (chrom == "chrY" & genderOfSamples[sam_no] == "F") {
        next
      } else if (chrom == "chrY" & genderOfSamples[sam_no] == "M") {
        initial_state <- 2
      } else {
        initial_state <- 3
      }
      start = left_borders[[l]]
      end = right_borders[[l]]
      for (k in 1:2) {
        output_of_plots <-  paste0(folder_name, sample_name)
        which_to_allow <- "NA"
        which_to_allow_ontarget <- "NA"
        if (k == 1) {
          which_to_allow = which(globalBed[,1] == chrom & as.numeric(globalBed[,2]) <= as.numeric(left_borders[[l]]) )
          which_to_allow_ontarget = which(bedFileFiltered[,1] == chrom & as.numeric(bedFileFiltered[,2]) <= as.numeric(left_borders[[l]]) )
        } else {
          which_to_allow = which(globalBed[,1] == chrom & as.numeric(globalBed[,2]) >= as.numeric(right_borders[[l]]) )
          which_to_allow_ontarget = which(bedFileFiltered[,1] == chrom & as.numeric(bedFileFiltered[,2]) >= as.numeric(right_borders[[l]]) )
        }
        if (length(which_to_allow) <= 1) {
          next
        }
        toyMatrixOfLikeliks = globalMatrixOfLikeliks[which_to_allow,]
        toybedFileFiltered = globalBed[which_to_allow,]
        found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, initial_state))
        
        # Due to inroduction of covariances intermediate probes we need to remap our variants back to the original bedFileFiltered
        if (covar & nrow(found_CNVs) > 0) {
          for (z in 1:nrow(found_CNVs)) {
            startOfFoundCNV = as.numeric(toybedFileFiltered[found_CNVs[z,2], 2])
            endOfFoundCNV = as.numeric(toybedFileFiltered[found_CNVs[z,3], 3])

            valuesInsideBed = which(bedFileFiltered[,1] == chrom & bedFileFiltered[,2] >= startOfFoundCNV & bedFileFiltered[,3] <= endOfFoundCNV)
            if (length(valuesInsideBed) > 0) {
              listForLikeliks = form_matrix_of_likeliks_one_sample_with_cov(1, ncol(coverage.normalised), sam_no, localSds[valuesInsideBed], coverage.normalised[valuesInsideBed,,drop=F], sqrt(cn_states / 2), covarianceTree, bedFileFiltered[valuesInsideBed,], threshold)
              if (!is.null(listForLikeliks)) {
                matrix_of_likeliks_with_covar_CNV <- listForLikeliks[[1]]
                bedFileFilteredWithArtificialProbesCNV = listForLikeliks[[2]]
                scoreToAdd = sum(matrix_of_likeliks_with_covar_CNV[,found_CNVs[z,4]] - matrix_of_likeliks_with_covar_CNV[,initial_state])
                found_CNVs[z,1] = as.numeric(found_CNVs[z,1]) + scoreToAdd
              }
            }
            
          }

        }

        toyCoverageGermline = coverage.normalised[which_to_allow_ontarget,sam_no]
        toySizesOfPointsFromLocalSds = sizesOfPointsFromLocalSds[which_to_allow]
        
        toyCoverageGermlineCohort = coverage.normalised[which_to_allow_ontarget,]

        
        if (nrow(found_CNVs) > 0) {
          alleleFrequency = rep(1 / ncol(coverage.normalised), nrow(found_CNVs))
          for (i in 1:nrow(found_CNVs)) {
            whichOnTarget = which(as.numeric(bedFileFiltered[which_to_allow_ontarget,2]) >= as.numeric(toybedFileFiltered[found_CNVs[i,2],2]) &
                                      as.numeric(bedFileFiltered[which_to_allow_ontarget,3]) <= as.numeric(toybedFileFiltered[found_CNVs[i,3],3])
                                      )
            cnState = cn_states[found_CNVs[i,4]]
            if (length(whichOnTarget) > 0) {
            mediansOfCoveragesInsideTheCohort <- apply(toyCoverageGermlineCohort[whichOnTarget,,drop=F], 2, median)
            if (cnState < 2) {
              alleleFrequency[i] = length(which(mediansOfCoveragesInsideTheCohort < (1 - (1 - sqrt(1/2)) / 2))) / ncol(coverage.normalised)
            }
            if (cnState > 2) {
              alleleFrequency[i] = length(which(mediansOfCoveragesInsideTheCohort > (1 + (sqrt(3/2) - 1) / 2))) / ncol(coverage.normalised)
            }
            } else {
              alleleFrequency[i] = -1.0
            }
          }
        }
        
        if (!chrom %in% c("chrX", "chrY")) {
          vectorOfZScoresLocaL = vectorOfZScores[which_to_allow]
          if (length(vectorOfZScoresLocaL) > 10)
            vectorWithNumberOfOutliers = c(vectorWithNumberOfOutliers, 
                                         length(which(vectorOfZScoresLocaL > qnorm(0.975) | vectorOfZScoresLocaL < qnorm(0.025))) / length(vectorOfZScoresLocaL))
           }
        
        ### IGV PLOTTING
        if (opt$visulizationIGV) {
          if(opt$debug) {
            print("START OF IGV PLOTTING")
          }
          if (sam_no_off) {
            toyCoverageGermline = c(coverage.normalised[,sam_no], coverage.normalised.off[,sam_no_off])[orderOfBed][which_to_allow]
          }
          
          outputFileNameCNVs <- paste0(folder_name, sample_name, "/", sample_name, "_cnvs.seg")
          outputFileNameDots <- paste0(folder_name, sample_name, "/", sample_name, "_cov.seg")
          reverseFunctionUsedToTransform = function(x, chrom) {return((2 * x ** 2))}
          outputSegmentsAndDotsFromListOfCNVs(toybedFileFiltered, found_CNVs, start, end, outputFileNameCNVs, 
                                              outputFileNameDots, sample_name, toyCoverageGermline, reverseFunctionUsedToTransform, cn_states)
          if(opt$debug) {
            print("END OF IGV PLOTTING")
          }
        }
        ### END OF IGV PLOTTING
        
        
        
        if (nrow(found_CNVs) > 0) {
          # UNCOMMENT FOR PLOTTING!!!
          
          cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyCoverageGermline, toybedFileFiltered, output_of_plots, chrom, cn_states, 
                                          toySizesOfPointsFromLocalSds,alleleFrequency, plottingOfPNGs)
          if (found_CNVs[1,1] != -1000) {
            found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
          }
          for (i in 1:nrow(found_CNVs)) {
            
            CNVnamesInside <- unlist(unique(toybedFileFiltered[found_CNVs[i,2]:found_CNVs[i,3],4]))
            if(opt$debug) {
              print(CNVnamesInside)
            }
            
            CNVentry = matrix(c(sample_name, chrom, toybedFileFiltered[found_CNVs[i,2],2], toybedFileFiltered[found_CNVs[i,3],3], 
                                paste(CNVnamesInside, collapse=", "),
                                found_CNVs[i,4] - 1, 
                                found_CNVs[i,5]),
                              nrow=1)
            if(opt$debug) {
              print(CNVentry)
            }
          }
        }
      }
      if (nrow(found_CNVs_total) > opt$maxNumGermCNVs & iterations != maxIteration) {
        break
      }
    }
    if (nrow(found_CNVs_total) < opt$maxNumGermCNVs) {
      numberOfCNVsIsSufficientlySmall = T
      break
    } else {
      if (iterations < maxIteration) {
        print("Sample had too many CNVs. We re-analyse it with stricter thresholds")
        threshold = threshold + 20
        minimum_length_of_CNV = minimum_length_of_CNV + 1
        unlink(paste0(folder_name, sample_name), recursive = T)
        dir.create(paste0(folder_name, sample_name))
      }
    }
  }
  
  ### FDR
  if (as.numeric(opt$fdrGermline) != 0) {
    positionsToExclude = c()
    for (z in 1:nrow(found_CNVs_total)) {
      if (as.numeric(found_CNVs_total[z,5] > 200)) {
        positionsToExclude = c(positionsToExclude, which(bedFileFiltered[,1] == found_CNVs_total[z,1] & as.numeric(bedFileFiltered[,2]) >= as.numeric(found_CNVs_total[z,2])
                                                                                                                                      & as.numeric(bedFileFiltered[,3]) <= as.numeric(found_CNVs_total[z,3])))
      }
    }
    matrix_of_likeliks_for_FDR = matrix_of_likeliks_for_FDR[-union( which(bedFileFiltered[,1] %in% c("chrX", "chrY")) ,  positionsToExclude),]
    numberOfIterationsForFDR = as.numeric(opt$fdrGermline)
    print(paste("Started to perform FDR permutations", Sys.time()))
    stepLikeChromosome = floor(nrow(matrix_of_likeliks_for_FDR) / 22)
    for (z in seq(from=1, to=nrow(matrix_of_likeliks_for_FDR) - stepLikeChromosome - 1, by=stepLikeChromosome)) {
      matrix_of_likeliks_for_FDR_part = matrix_of_likeliks_for_FDR[z:(z + stepLikeChromosome),]
      detectedFalseCNVs <- foreach(i=1:numberOfIterationsForFDR, .combine="rbind") %dopar% {
        shuffledMatrixOfLikelis = matrix_of_likeliks_for_FDR_part[sample(1:nrow(matrix_of_likeliks_for_FDR_part)),1:(main_initial_state + 2)]
        detectedCnvs <- find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, main_initial_state, shuffledMatrixOfLikelis, main_initial_state)
        detectedCnvs
      }
    }
    print(paste("Fnished to perform FDR permutations", Sys.time()))
    detectedDeletions <-  detectedFalseCNVs[which(detectedFalseCNVs[,4] < main_initial_state),1:4,drop=F]
    detectedDuplications <-  detectedFalseCNVs[which(detectedFalseCNVs[,4] > main_initial_state),1:4,drop=F]
    
    
    if (nrow(detectedDeletions) > 0) {
      thresholdsDel = sort(-1 * unique(detectedDeletions[,1]))
      fdrThreshold = 0.05
      currentThresholdDel = thresholdsDel[i]
      for (i in 1:length(thresholdsDel)) {
        currentThresholdDel = thresholdsDel[i]
        FDR = length(which(-1 * detectedDeletions[,1] > currentThresholdDel)) / (length(which(!found_CNVs_total[,1] %in% c("chrX", "chrY") & as.numeric(found_CNVs_total[,4]) < 2 & as.numeric(found_CNVs_total[,5]) > currentThresholdDel)) )
        if (FDR < fdrThreshold) {
          break
        }
      }
      currentThresholdDel = currentThresholdDel + 1
    } else {
      currentThresholdDel = threshold
    }
    
    if (nrow(detectedDuplications) > 0) {
      thresholdsDup = sort(-1 * unique(detectedDuplications[,1]))
      fdrThreshold = 0.05
      currentThresholdDup = thresholdsDup[i]
      for (i in 1:length(thresholdsDup)) {
        currentThresholdDup = thresholdsDup[i]
        FDR = length(which(-1 * detectedDuplications[,1] > currentThresholdDel)) / (length(which(!found_CNVs_total[,1] %in% c("chrX", "chrY") & as.numeric(found_CNVs_total[,4]) > 2 & as.numeric(found_CNVs_total[,5]) > currentThresholdDup)) )
        if (FDR < fdrThreshold) {
          break
        }
      }
      currentThresholdDup = currentThresholdDup + 1
    } else {
      currentThresholdDup = threshold
    }
    
    columnOfFilter = matrix("NOT PASS", ncol=1, nrow=nrow(found_CNVs_total))
    for (i in 1:nrow(found_CNVs_total)) {
      if (as.numeric(found_CNVs_total[i,4]) < 2) {
        if (as.numeric(found_CNVs_total[i,5]) > currentThresholdDel) columnOfFilter[i] = "PASS"
      } else {
        if (as.numeric(found_CNVs_total[i,5]) > currentThresholdDup) columnOfFilter[i] = "PASS"
      }
    }
    found_CNVs_total = cbind(found_CNVs_total, columnOfFilter)
    colnames(found_CNVs_total)[ncol(found_CNVs_total)] = "FDR_filter"
  }
  
  finalPValue = 1.0
  fileToOut <- paste0(folder_name, sample_name, paste0("/", sample_name, "_cnvs.tsv"))
  fileConn<-file(fileToOut)
  writeLines(c(paste("##number of iterations:", iterations, ", gender of sample:", genderOfSamples[sam_no], ", was it outlier after clustering?", outliersByClustering[sam_no], collapse = " "), 
               paste("##fraction of outliers:", round(median(vectorWithNumberOfOutliers), digits=3), collapse = " ")), fileConn)
  close(fileConn)
  found_CNVs_total[,7] = (format(as.numeric(found_CNVs_total[,7]), nsmall=3))
  found_CNVs_total[,8] = (format(as.numeric(found_CNVs_total[,8]), nsmall=3))
  write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)
}

stopCluster(cl)
