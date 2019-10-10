library(party)



cn_states <- 0:8
if (opt$mosaicism) {
  cn_states_mosaicism <- seq(from=1.1, to=2.9, by=0.05)
  
  diffs <- sapply(cn_states_mosaicism, function(i) {min(abs(cn_states - i))})
  cn_states_mosaicism = cn_states_mosaicism[which(diffs > 0.09 & diffs < 0.91)]
  cn_states_mosaicism = unique(c(cn_states_mosaicism, cn_states))
}



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
if (ncol(coverage.normalised) > 25) {
  covar = T
} else {
  covar = F
}

print(paste("We start to estimate covariances between neighboring regions in germline data - may take some time", Sys.time()))
setwd(opt$out)
covarianceTree = NULL
if (covar == T) {
  covarianceTree = returnTreeForCorrelation(coverage.normalised, 
                                            sdsOfGermlineSamples, sdsOfProbes, 
                                            bedFileFiltered)
}
if (is.null(covarianceTree)) {
  covar = F
}
setwd(opt$folderWithScript)
print(paste("Tree of covariances (using 2 predictors - sum of regions' lengths and log2 of distance between regions) plotted in", opt$out, Sys.time()))

positionsInPolymorphic = c()
if (!is.null(polymorphicRegions)) {
  for (chrom in unique(bedFileFiltered[,1])) {
    chromInBed = which(bedFileFiltered[,1] == chrom)
    polymorphicInsideChrom = polymorphicRegions[which(polymorphicRegions[,1] == chrom),]
    for (i in 1:nrow(polymorphicInsideChrom)) {
      whichInsideVariant <- which(as.numeric(bedFileFiltered[chromInBed,2]) >= as.numeric(polymorphicInsideChrom[i,2]) - 500 & as.numeric(bedFileFiltered[chromInBed,3]) <= as.numeric(polymorphicInsideChrom[i,3]) + 500)
      positionsInPolymorphic = c(positionsInPolymorphic, chromInBed[whichInsideVariant])
    }
  }
}

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
  sdsForOutput = localSds
  
  if (sam_no_off) {
    localSdsOff = sdsOfProbesOff * sdsOfGermlineSamplesOff[sam_no_off]
    localSdsOff[which(localSdsOff == 0)] = median(localSdsOff)
  }
  
  
  if (!dir.exists(paste0(folder_name, sample_name))) {
    dir.create(paste0(folder_name, sample_name))
  } else {
    if (is.null(opt$reanalyseCohort)) next
  }
  setwd(paste0(folder_name, sample_name))
  
  
  dict_to_output = c()
  
  
  matrix_of_likeliks_for_FDR <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised), sam_no, localSds, coverage.normalised, sqrt(cn_states / 2))
  if (opt$mosaicism) {
    matrix_of_likeliks_for_FDR_mosaic <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised), sam_no, localSds, coverage.normalised, sqrt(cn_states_mosaicism / 2))
  }
  if (!is.null(polymorphicRegions)) {
    matrix_of_likeliks_for_FDR[positionsInPolymorphic,] = 0
    matrix_of_likeliks_for_FDR[positionsInPolymorphic,main_initial_state] = 0
    if (opt$mosaicism) {
      matrix_of_likeliks_for_FDR_mosaic[positionsInPolymorphic,] = 0
      matrix_of_likeliks_for_FDR_mosaic[positionsInPolymorphic,main_initial_state] = 0
    }
  }
  if (opt$mosaicism) {
    fineForMosaicism = 0.05
    matrix_of_likeliks_for_FDR_mosaic[,which(cn_states_mosaicism %% 1 != 0)] = matrix_of_likeliks_for_FDR_mosaic[,which(cn_states_mosaicism %% 1 != 0)] + fineForMosaicism
  }
  matrix_of_likeliks <- matrix_of_likeliks_for_FDR
  
  
  
  if (sam_no_off) {
    matrix_of_likeliks_off = form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised.off), sam_no_off, localSdsOff, coverage.normalised.off, sqrt(cn_states / 2))
    if (opt$mosaicism) {
      matrix_of_likeliks_off_mosaic = form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised.off), sam_no_off, localSdsOff, coverage.normalised.off, sqrt(cn_states_mosaicism / 2))
    }
    globalBed = rbind(bedFileFiltered, bedFileFilteredOfftarget)
    orderOfBed = order(globalBed[,1], as.numeric(globalBed[,2]))
    globalBed = globalBed[orderOfBed,]
    globalMatrixOfLikeliks = rbind(matrix_of_likeliks, matrix_of_likeliks_off)
    globalMatrixOfLikeliks = globalMatrixOfLikeliks[orderOfBed,]
    if (opt$mosaicism) {
      globalMatrixOfLikeliksMosaic = rbind(matrix_of_likeliks_for_FDR_mosaic, matrix_of_likeliks_off_mosaic)[orderOfBed,]
    }
    sdsForOutput = c(localSds, localSdsOff)[orderOfBed]
  } else {
    globalBed = bedFileFiltered
    globalMatrixOfLikeliks = matrix_of_likeliks
    if (opt$mosaicism) {
      globalMatrixOfLikeliksMosaic = matrix_of_likeliks_for_FDR_mosaic
    }
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
        armFinalized = F
        local_cn_states = cn_states
        presenceOfMosaicVariants = F
        while(!armFinalized) {
          
          if (opt$mosaicism == F | chrom %in% c("chrX", "chrY")) armFinalized = T
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
            armFinalized = T
            next
          }
          polymorphicFromThisArm = c()
          if (length(polymorphicRegions) > 0 & opt$mosaicism) {
            polymorphicFromThisArm = intersect(polymorphicRegions, which_to_allow_ontarget)
          }
          toyMatrixOfLikeliks = globalMatrixOfLikeliks[which_to_allow,]
          toybedFileFiltered = globalBed[which_to_allow,]
          if (presenceOfMosaicVariants == F) {
            found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, initial_state))
          } else {
            local_cn_states = cn_states_mosaicism
            toyMatrixOfLikeliksMosaic = globalMatrixOfLikeliksMosaic[which_to_allow,]
            found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliksMosaic, initial_state))
            found_CNVs_mosaic_too_short = which(local_cn_states[found_CNVs[,4]] %% 1 != 0 & (found_CNVs[,3]- found_CNVs[,2]) < max(10, 3 * opt$lengthG))
            if (length(found_CNVs_mosaic_too_short) > 0) {
              found_CNVs = found_CNVs[-found_CNVs_mosaic_too_short,,drop=F]
            }
            armFinalized = T
          }
          found_CNVs_recall = NULL
          if (opt$superRecall < opt$scoreG) {
            if (nrow(found_CNVs) > 0) {
              for (z in 1:nrow(found_CNVs)) {
                startOfCNV = found_CNVs[z,2]
                endOfCNV = found_CNVs[z,3]
                toyMatrixOfLikeliks[max(1, startOfCNV - 1):min(endOfCNV + 1, nrow(toyMatrixOfLikeliks)),] = 0
              }
            }
            found_CNVs_recall <- as.matrix(find_all_CNVs(0, opt$superRecall, 0, initial_state, toyMatrixOfLikeliks, initial_state))
          }
          if (!is.null(found_CNVs_recall)) {
            found_CNVs = rbind(found_CNVs, found_CNVs_recall)
          }
          
          
          # Due to inroduction of covariances intermediate probes we need to remap our variants back to the original bedFileFiltered
          if (covar & nrow(found_CNVs) > 0) {
            for (z in 1:nrow(found_CNVs)) {
              startOfFoundCNV = as.numeric(toybedFileFiltered[found_CNVs[z,2], 2])
              endOfFoundCNV = as.numeric(toybedFileFiltered[found_CNVs[z,3], 3])
              
              valuesInsideBed = which(bedFileFiltered[,1] == chrom & bedFileFiltered[,2] >= startOfFoundCNV & bedFileFiltered[,3] <= endOfFoundCNV)
              if (length(valuesInsideBed) > 1) {
                listForLikeliks = form_matrix_of_likeliks_one_sample_with_cov(1, ncol(coverage.normalised), sam_no, localSds[valuesInsideBed], coverage.normalised[valuesInsideBed,,drop=F], sqrt(local_cn_states / 2), covarianceTree, bedFileFiltered[valuesInsideBed,], threshold)
                if (!is.null(listForLikeliks)) {
                  matrix_of_likeliks_with_covar_CNV <- listForLikeliks[[1]]
                  bedFileFilteredWithArtificialProbesCNV = listForLikeliks[[2]]
                  scoreToAdd = sum(matrix_of_likeliks_with_covar_CNV[,found_CNVs[z,4]] - matrix_of_likeliks_with_covar_CNV[,initial_state])
                  found_CNVs[z,1] = as.numeric(found_CNVs[z,1]) + scoreToAdd
                }
              }
              
            }
            
          }
          
          toySizesOfPointsFromLocalSds = sizesOfPointsFromLocalSds[which_to_allow]
          toyCoverageGermline = coverage.normalised[which_to_allow_ontarget,sam_no]
          toyCoverageGermlineCohort = coverage.normalised[which_to_allow_ontarget,]
          
          ### CHECKING FOR MOSAICISM!
          if (length(which_to_allow_ontarget) > max(10, opt$lengthG * 2)) {
            positionsToRemove = c()
            if (nrow(found_CNVs) > 0) {
              for (z in 1:nrow(found_CNVs)) {
                positionsToRemove = c(positionsToRemove, found_CNVs[z,2]:found_CNVs[z,3])
              }
            }
            if (length(polymorphicFromThisArm) > 0) {
              positionsToRemove = c(positionsToRemove, polymorphicFromThisArm)
            }
            toyCoverageGermlineWithoutNonMosaicCNVs = toyCoverageGermline
            if (length(positionsToRemove) > 0) {
              toyCoverageGermlineWithoutNonMosaicCNVs = toyCoverageGermlineWithoutNonMosaicCNVs[-positionsToRemove]
            }
            if (length(toyCoverageGermlineWithoutNonMosaicCNVs) <= max(10, opt$lengthG * 2)) {print("Less than required number of regions for mosaic CNVs detection");armFinalized = T} else {
              rollingThrowCoverage = runmed(toyCoverageGermlineWithoutNonMosaicCNVs, max(10, opt$lengthG * 2))
              standDevOfRolling = Qn(rollingThrowCoverage)
              locationOfRolling = median(rollingThrowCoverage)
              if (max(abs(rollingThrowCoverage - 1)) > max(3 * standDevOfRolling, 0.05) & !chrom %in% c("chrX", "chrY")) {
                print(paste("Potential mosaicism at chromosome", chrom, ", arm: ", k, "(this is normal if you have not filtered polymorphic regions before calling)"))
                if (opt$mosaicism) {
                  presenceOfMosaicVariants = T
                  if (!armFinalized)
                    next
                }
              } else {
                armFinalized = T
              }
            }
          } else {
            armFinalized = T
          }
          ### CHECKING FOR MOSAICISM OVER!
          
          if (nrow(found_CNVs) > 0) {
            alleleFrequency = rep(1 / ncol(coverage.normalised), nrow(found_CNVs))
            for (i in 1:nrow(found_CNVs)) {
              whichOnTarget = which(as.numeric(bedFileFiltered[which_to_allow_ontarget,2]) >= as.numeric(toybedFileFiltered[found_CNVs[i,2],2]) &
                                      as.numeric(bedFileFiltered[which_to_allow_ontarget,3]) <= as.numeric(toybedFileFiltered[found_CNVs[i,3],3])
              )
              cnState = local_cn_states[found_CNVs[i,4]]
              if (chrom %in% c("chrX", "chrY")) {
                allowedSamples <- which(genderOfSamples == genderOfSamples[sam_no])
              } else {
                allowedSamples = 1:ncol(toyCoverageGermlineCohort)
              }
              if (length(whichOnTarget) > 0) {
                mediansOfCoveragesInsideTheCohort <- apply(toyCoverageGermlineCohort[whichOnTarget,allowedSamples,drop=F], 2, median)
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
                                                outputFileNameDots, sample_name, toyCoverageGermline, reverseFunctionUsedToTransform, local_cn_states, sdsForOutput[which_to_allow])
            if(opt$debug) {
              print("END OF IGV PLOTTING")
            }
          }
          ### END OF IGV PLOTTING
          
          
          
          if (nrow(found_CNVs) > 0) {
            # UNCOMMENT FOR PLOTTING!!!
            
            cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyCoverageGermline, toybedFileFiltered, output_of_plots, chrom, local_cn_states, 
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
      }
      if (length(which(as.numeric(found_CNVs_total[,5]) > threshold)) > opt$maxNumGermCNVs & iterations != maxIteration) {
        break
      }
    }
    if (length(which(as.numeric(found_CNVs_total[,5]) > threshold))  < opt$maxNumGermCNVs) {
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
  
  pvaluesForCNVs <- rep(NA, nrow(found_CNVs_total))
  if (nrow(found_CNVs_total) > 0) {
    for (i in 1:nrow(found_CNVs_total)) {
      coordsInBedOn = which(bedFileFiltered[,1] == found_CNVs_total[i,1] & as.numeric(bedFileFiltered[,2]) >= as.numeric(found_CNVs_total[i,2]) & as.numeric(bedFileFiltered[,3]) <= as.numeric(found_CNVs_total[i,3]))
      if (length(coordsInBedOn) > 0) {
        valueOfSample = median(coverage.normalised[coordsInBedOn,sam_no])
        if (!chrom %in% c("chrX","chrY")) {
          valueOfOthers = apply(coverage.normalised[coordsInBedOn,-sam_no,drop=F],2,median)
        } else {
          valueOfOthers = apply(coverage.normalised[coordsInBedOn,which(genderOfSamples == genderOfSamples[sam_no]),drop=F],2,median)
        }
        sdOfOthers = Qn(valueOfOthers)
        pval = 2 * (pnorm(-abs(valueOfSample - median(valueOfOthers)) / sdOfOthers))
        pvaluesForCNVs[i] = pval
      } else {
        if (sam_no_off > 0) {
          coordsInBedOff = which(bedFileFilteredOfftarget[,1] == found_CNVs_total[i,1] & as.numeric(bedFileFilteredOfftarget[,2]) >= as.numeric(found_CNVs_total[i,2]) & as.numeric(bedFileFilteredOfftarget[,3]) <= as.numeric(found_CNVs_total[i,3]))
          if (length(coordsInBedOff) > 0) {
            valueOfSample = median(coverage.normalised.off[coordsInBedOff,sam_no_off])
            if (!chrom %in% c("chrX","chrY")) {
              valueOfOthers = apply(coverage.normalised.off[coordsInBedOff,-sam_no,drop=F],2,median)
            } else {
              valueOfOthers = apply(coverage.normalised.off[coordsInBedOff,which(genderOfSamples == genderOfSamples[sam_no]),drop=F],2,median)
            }
          }
        }
        sdOfOthers = Qn(valueOfOthers)
        pval = 2 * (pnorm(-abs(valueOfSample - median(valueOfOthers)) / sdOfOthers))
        pvaluesForCNVs[i] = pval
      }
    }
    pvaluesForCNVs = p.adjust(pvaluesForCNVs, method="fdr")
    found_CNVs_total = cbind(found_CNVs_total, format(round(pvaluesForCNVs, 5), scientific = F))
    colnames(found_CNVs_total)[ncol(found_CNVs_total)] = "qvalue"
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
  
  if (opt$mosaicism) {
    # BLOCK WITH CLONALITY
  }
  
  finalPValue = 1.0
  fileToOut <- paste0(folder_name, sample_name, paste0("/", sample_name, "_cnvs.tsv"))
  fileConn<-file(fileToOut)
  writeLines(c(
    paste0("##ANALYSISTYPE=CLINCNV_GERMLINE_SINGLE"), 
    paste0("##", clincnvVersion), 
    paste("##Analysis finished on:", Sys.time()),
    paste("##gender of sample:", genderOfSamples[sam_no], collapse = " "),
    paste("##number of iterations:", iterations, collapse = " "), 
    paste("##quality used at final iteration:", threshold, collapse = " "), 
    paste("##was it outlier after clustering:", outliersByClustering[sam_no], collapse = " "),
    paste("##fraction of outliers:", round(median(vectorWithNumberOfOutliers), digits=3), collapse = " ")), fileConn)
  close(fileConn)
  found_CNVs_total[,7] = (format(as.numeric(found_CNVs_total[,7]), nsmall=3))
  found_CNVs_total[,8] = (format(as.numeric(found_CNVs_total[,8]), nsmall=3))
  write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)
}


