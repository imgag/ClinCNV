
somaticCalling <- function(matrixOfLogFold) {
  for (sam_no in 1:ncol(matrixOfLogFold)) {

    clusterExport(cl,c('maxSubArraySum', 'fillInPList', 'likelihoodOfSNV','return_likelik', 'vect_of_norm_likeliks', 'vect_of_t_likeliks'))
    sample_name <- colnames(matrixOfLogFold)[sam_no]
    overdispersionNormal = NULL
    sampleInOfftarget = F
    local_datasetOfPuritiesCopies = datasetOfPuritiesCopies
    local_datasetOfPuritiesCopiesForFinal = datasetOfPuritiesCopiesForFinalIteration
    
    germline_sample_no = which(colnames(tmpNormal) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2])
    germline_sample_name = colnames(tmpNormal)[germline_sample_no]
    tumor_sample_no = which(colnames(tumor) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1])
    if (!frameworkOff == "ontarget")
      if (strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1] %in% colnames(tmpTumorOff)) {
        tumor_sample_no_off = which(colnames(tmpTumorOff) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1])
        tumor_sample_name = strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
      }
    if (germline_sample_no == "F") next
    
    if (!is.null(opt$normalSample) & !is.null(opt$tumorSample)) {
      if (!sample_name == paste(opt$tumorSample, opt$normalSample, sep="-")) {
        next
      }
    }

    # To speed up reiteration, we do not want match between BAF file and bed file a lot of times
    if (frameworkDataTypes == "covdepthBAF") {
      closestBedRegions <- c()
      vectorsWithRegionCoordsFilled = F
    }
    
    print(paste("We are working on sample name:", sample_name))
    print(Sys.time())
    
    
    
    if (frameworkOff == "offtarget") {
      if (sample_name %in% colnames(matrixOfLogFoldOff)) {
        sam_no_off = which(colnames(matrixOfLogFoldOff) == sample_name)
        sampleInOfftarget = T
      }
    }
    
    
    matrixOfClonality = matrix(0, nrow=1, ncol=1)
    
    
    
    if (!dir.exists(paste0(folder_name, sample_name)) | (!is.null(opt$reanalyseCohort))) {
      
      dir.create(paste0(folder_name, sample_name))
      setwd(paste0(folder_name, sample_name))
      
      
      finalIteration = F
      shiftToTry = opt$shiftToTry
      while(T) {
        #### CORRECTION - IF THE SAMPLE HAS TOO MANY CNAS, WE EXPECT SOME SHIFT THERE
        if (!finalIteration) {
          if (frameworkDataTypes == "covdepthBAF" & germline_sample_name %in% normalNames) {
            sampleName2 <- strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
            tumorNames = sapply(1:length(allowedChromsBaf), function(i) {strsplit(names(allowedChromsBaf)[i], split="-")[[1]][1]})
            position <- which(tumorNames == sampleName2)
            if (length(position) == 1) {
              allowedChromsBafSample <- allowedChromsBaf[[position]]
              if (!sampleInOfftarget) {
                shiftOfCoverage = find_baseline_level(allowedChromsBafSample, matrixOfLogFold[,sam_no], bedFileForCluster)  
                print(paste("Potential shift lines of diploid states", paste(round(shiftOfCoverage, digits=4), collapse="; "), " - if there are more than 1 different shift, you may specify which one you want to choose by putting --shiftToTry flag"))
                if (length(shiftOfCoverage) >= shiftToTry) {
                  shiftOfCoverage = shiftOfCoverage[shiftToTry]
                }
              } else {
                sam_no_off = which(colnames(matrixOfLogFoldOff) == sample_name)
                shiftOfCoverage = find_baseline_level(allowedChromsBafSample, matrixOfLogFold[,sam_no], bedFileForCluster, matrixOfLogFoldOff[,sam_no_off], bedFileForClusterOff)   
                print(paste("Potential shift lines of diploid states", paste(round(shiftOfCoverage, digits=4), collapse="; "), " - if there are more than 1 different shift, you may specify which one you want to choose by putting --shiftToTry flag"))
                
                if (length(shiftOfCoverage) >= shiftToTry) {
                  shiftOfCoverage = shiftOfCoverage[shiftToTry]
                }
              }
              if (length(shiftOfCoverage) == 1) {
                matrixOfLogFold[,sam_no] = matrixOfLogFold[,sam_no] - shiftOfCoverage
                if (sampleInOfftarget)
                  matrixOfLogFoldOff[,sam_no_off] = matrixOfLogFoldOff[,sam_no_off] - shiftOfCoverage
              } else {
                print("Several shifts of coverage specified! We do not apply any shift of coverage.")
              }
            }
          }
        }
        # CLEAN FOLDER IN THE BEGINNING OF EACH ITERATION
        if (!finalIteration)
          do.call(file.remove, list(list.files(paste0(folder_name, sample_name), full.names = TRUE)))
        
        if (finalIteration ) {
          indices_to_remove_by_purity <- which(!(local_datasetOfPuritiesCopiesForFinal[,6] %in% c(0, clonalBestPurities)))
          if (length(indices_to_remove_by_purity) > 0)
            local_datasetOfPuritiesCopies = local_datasetOfPuritiesCopiesForFinal[-indices_to_remove_by_purity,]
          
          
          potential_second_copy_changes <- c(unique(dist(unique(local_datasetOfPuritiesCopies[,6]))))
          #potential_second_copy_changes <- c(unique(local_datasetOfPuritiesCopies[,6]))
          potential_second_copy_changes <- c(setdiff(potential_second_copy_changes, max(potential_second_copy_changes)), 0)
          if (length(potential_second_copy_changes) > 0) {
            indices_to_remove_by_purity <- which(!(local_datasetOfPuritiesCopies[,9] %in% c(0, potential_second_copy_changes)))
            if (length(indices_to_remove_by_purity) > 0)
              local_datasetOfPuritiesCopies = local_datasetOfPuritiesCopies[-indices_to_remove_by_purity,]
          }
        }
        
        local_purities <- local_datasetOfPuritiesCopies[,6]
        local_copy_numbers_used_major <- local_datasetOfPuritiesCopies[,2]
        local_copy_numbers_used_minor <- local_datasetOfPuritiesCopies[,3]
        local_cn_states <- local_datasetOfPuritiesCopies[,1]
        local_majorBAF <- local_datasetOfPuritiesCopies[,4]
        local_minorBAF <- local_datasetOfPuritiesCopies[,5]
        local_copy_numbers_used_major_second <- local_datasetOfPuritiesCopies[,7]
        local_copy_numbers_used_minor_second <- local_datasetOfPuritiesCopies[,8]
        local_purities_second <- local_datasetOfPuritiesCopies[,9]
        
        # PART FOR MATRIX OF CLONALITY (ONLY 2 CLONES)
        uniqueLocalPurities = unique(local_purities)
        zeroPurity <- which(uniqueLocalPurities == 0)
        if (length(zeroPurity) > 0)
          uniqueLocalPurities = sort(uniqueLocalPurities[-zeroPurity])
        uniqueLocalPuritiesSecond = unique(local_purities_second)
        likeliksFoundCNVsVsPuritiesGlobal = matrix(nrow=0, ncol=length(uniqueLocalPurities) * length(uniqueLocalPuritiesSecond))
        
        
        pvalsForQC <- c()
        threshold = opt$scoreS
        minimum_length_of_CNV = opt$lengthS
        
        #if (!finalIteration) {
        #  threshold = opt$scoreS + 100
        #}
        
        price_per_tile = 0.025
        initial_state <- 1
        
        
        localSds = sdsOfProbes * (sdsOfSomaticSamples[sam_no])
        if (frameworkOff == "offtarget") {
          if (sampleInOfftarget) {
            localSdsOff = sdsOfProbesOff * (sdsOfSomaticSamplesOff[sam_no_off])
          }
        }
        
        
        
        dict_to_output = c()
        
        
        
        
        
        
        
        
        matrixOfLogFoldCorrectedSmall = matrixOfLogFold
        matrixOfLogFoldCorrectedSmall[which(matrixOfLogFoldCorrectedSmall < log2(min(local_cn_states) / 2))] = log2(min(local_cn_states) / 2)
        matrixOfLogFoldCorrectedSmall[which(matrixOfLogFoldCorrectedSmall > log2(max(local_cn_states) / 2))] = log2(max(local_cn_states) / 2)
        
        
        if (genderOfSamples[germline_sample_no] == "F") {
          whichAllowedOnTagret = which(bedFileForCluster[,1] != "chrY")
        } else {
          whichAllowedOnTagret = 1:nrow(bedFileForCluster)
        }
        if (sampleInOfftarget) {
          if (genderOfSamples[germline_sample_no] == "F") {
            whichAllowedOffTagret = which(bedFileForClusterOff[,1] != "chrY")
          } else {
            whichAllowedOffTagret = 1:nrow(bedFileForClusterOff)
          }
          blocked_states_global = determine_potential_states(matrixOfLogFoldCorrectedSmall[whichAllowedOnTagret,sam_no], local_cn_states, matrixOfLogFoldOff[whichAllowedOffTagret,sam_no_off])
        } else {
          blocked_states_global = determine_potential_states(matrixOfLogFoldCorrectedSmall[whichAllowedOnTagret,sam_no], local_cn_states)
        }
        
        if (length(blocked_states_global) > 0) {
          local_purities <- local_purities[-blocked_states_global]
          local_copy_numbers_used_major <- local_copy_numbers_used_major[-blocked_states_global]
          local_copy_numbers_used_minor <- local_copy_numbers_used_minor[-blocked_states_global]
          local_cn_states <- local_cn_states[-blocked_states_global]
          local_majorBAF <- local_majorBAF[-blocked_states_global]
          local_minorBAF <- local_minorBAF[-blocked_states_global]
          local_copy_numbers_used_major_second <- local_copy_numbers_used_major_second[-blocked_states_global]
          local_copy_numbers_used_minor_second <- local_copy_numbers_used_minor_second[-blocked_states_global]
          local_purities_second <- local_purities_second[-blocked_states_global]
        }
        
        #matrix_of_likeliks <- form_matrix_of_likeliks_one_sample(1,  matrixOfLogFoldCorrectedSmall[,sam_no], localSds, log2(local_cn_states/2))
        
        
        
        #if (genderOfSamples[germline_sample_no] == "M") {
        #  if (length(which(bedFileForCluster[,1] %in% c("chrX","chrY"))) > 0)
        #    matrix_of_likeliks[which(bedFileForCluster[,1] %in% c("chrX","chrY")),] = form_matrix_of_likeliks_one_sample(
        #      1, matrixOfLogFoldCorrectedSmall[which(bedFileForCluster[,1] %in% c("chrX","chrY")),sam_no], 
        #      localSds[which(bedFileForCluster[,1] %in% c("chrX","chrY"))], log2((1 - local_purities) + local_purities * local_copy_numbers_used_major))
        #}
        
        
        ### ADD LIKELIHOODS
        bAlleleFreqsTumor = NULL
        bAlleleFreqsNormal = NULL
        bafDeviationsForComparison = NULL
        matrixOfBAFLikeliks = NULL
        if (frameworkDataTypes == "covdepthBAF" & germline_sample_name %in% normalNames) {
          print("Started BAF calculation")
          print(Sys.time())
          
          numberOfAssignedPositions = 0
          sampleName2 <- strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
          tumorNames = sapply(1:length(allowedChromsBaf), function(i) {strsplit(names(allowedChromsBaf)[i], split="-")[[1]][1]})
          position <- which(tumorNames == sampleName2)
          if (length(position) == 1) {
            bAlleleFreqsTumor <- bAlleleFreqsAllSamples[[position]][[ strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1] ]]
            bAlleleFreqsNormal <- bAlleleFreqsAllSamples[[position]][[ strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2] ]]
            degreeOfRoughness = round(quantile(bAlleleFreqsTumor[,6], 0.5))
            if (finalIteration) {
              degreeOfRoughness = round(max(bAlleleFreqsTumor[,6]))
            }
            if (genderOfSamples[germline_sample_no] == "M") {
              bAlleleFreqsTumor = bAlleleFreqsTumor[which(!bAlleleFreqsTumor[,1] %in% c("chrX", "chrY")),]
              bAlleleFreqsNormal = bAlleleFreqsNormal[which(!bAlleleFreqsNormal[,1] %in% c("chrX", "chrY")),]
            }
            overdispersionNormal = overdispersionsNormal[[position]]
            overdispersionTumor = overdispersionsTumor[[position]]
            pvalueShift = pvaluesShifts[[position]]
            # calculate median correction factor
            allowedChromosomesAutosomesOnlySNV = c()
            allowedChromsBafSample <- allowedChromsBaf[[position]]
            for (allowedArm in allowedChromsBafSample) {
              splittedValue <- strsplit(allowedArm, "-")
              allowedChromosomesAutosomesOnlySNV = c(allowedChromosomesAutosomesOnlySNV, which(bAlleleFreqsTumor[,1] == splittedValue[[1]][1] & 
                                                                                                 as.numeric(bAlleleFreqsTumor[,2]) >= as.numeric(splittedValue[[1]][2]) &
                                                                                                 as.numeric(bAlleleFreqsTumor[,3]) <= as.numeric(splittedValue[[1]][3])
              )
              )
            }
            
            multiplierOfSNVsDueToMapping <- median(as.numeric(bAlleleFreqsTumor[allowedChromosomesAutosomesOnlySNV,5]))
            if (abs(multiplierOfSNVsDueToMapping - median(as.numeric(bAlleleFreqsNormal[,5]))) > 0.025) {
              multiplierOfSNVsDueToMapping = median(as.numeric(bAlleleFreqsNormal[,5]))
            }
            print("Multiplier of allele balance of a particular sample")
            print(multiplierOfSNVsDueToMapping)
            print("Multiplier of allele balance of its normal pair")
            print(median(as.numeric(bAlleleFreqsNormal[,5])))
            
            bafDeviationsForComparison = rep(0, length(allowedChromosomesAutosomesOnlySNV))
            counter = 1
            for (bafFromAllowedChr in allowedChromosomesAutosomesOnlySNV) {
              diff = (as.numeric(bAlleleFreqsTumor[bafFromAllowedChr,5]) - multiplierOfSNVsDueToMapping) * as.numeric(bAlleleFreqsTumor[bafFromAllowedChr,6]) 
              sds = sqrt(overdispersionTumor[bafFromAllowedChr] * as.numeric(bAlleleFreqsTumor[bafFromAllowedChr,6]) * multiplierOfSNVsDueToMapping * ( 1 - multiplierOfSNVsDueToMapping))
              bafDeviationsForComparison[counter] = abs(diff / sds)
              counter = counter + 1
            }
            
            
            if (length(closestBedRegions) == 0) closestBedRegions = rep(0, nrow(bAlleleFreqsTumor))
            for (i in 1:nrow(bAlleleFreqsTumor)) {
              # To avoid computationally expensive steps on the start of estimation
              if (vectorsWithRegionCoordsFilled) {
                closestBedRegion = closestBedRegions[i]
              } else {
                closestBedRegion <- which(bedFileForCluster[,1] == bAlleleFreqsTumor[i,1] & as.numeric(bedFileForCluster[,2]) - 50 <= as.numeric(bAlleleFreqsTumor[i,2])
                                          & as.numeric(bedFileForCluster[,3]) + 50 >= as.numeric(bAlleleFreqsTumor[i,3]))
                if (length(closestBedRegion) >= 1) {
                  closestBedRegion = closestBedRegion[1]
                  closestBedRegions[i] = closestBedRegion
                } else {
                  # if there is no matching position we put 0
                  closestBedRegions[i] = 0
                  closestBedRegion = 0
                }
              }
            }
            vectorsWithRegionCoordsFilled = T
            #bAlleleFreqsTumor = bAlleleFreqsTumor[which(closestBedRegions != 0),]
            #closestBedRegions = closestBedRegions[which(closestBedRegions != 0)]
            noBAF = F
            if (!finalIteration) {
              
              pvalues <- rep(0, nrow(bAlleleFreqsTumor))
              for (l in 1:nrow(bAlleleFreqsTumor)) {
                refAleleTum <- (round(as.numeric(bAlleleFreqsTumor[l,5]) * as.numeric(bAlleleFreqsTumor[l,6])))
                altAleleTum <- as.numeric(bAlleleFreqsTumor[l,6]) - refAleleTum
                refAleleNorm<- (round(as.numeric(bAlleleFreqsNormal[l,5]) * as.numeric(bAlleleFreqsNormal[l,6])))
                altAleleNorm <- as.numeric(bAlleleFreqsNormal[l,6]) - refAleleNorm
                pvalues[l] <- passPropTest(refAleleTum, refAleleNorm, altAleleTum, altAleleNorm)
              }
              pvaluesSmall = filter(pvalues < 0.01, rep(1 / 20, 20), sides = 2)
              pvaluesMed = runmed(pvalues, 5)
              
              # PRESEGMENTATION
              unique_local_cn_states = unique(local_cn_states)
              matrix_of_likeliks_for_presegmentation <- form_matrix_of_likeliks_one_sample(1,  matrixOfLogFoldCorrectedSmall[which(!bedFileForCluster[,1] %in% c("chrX","chrY")),sam_no], localSds, log2(unique_local_cn_states/2))
              presegmentation <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, matrix_of_likeliks_for_presegmentation, 1))
              notAllowedChromosomesAutosomesCNVs = c()
              if (nrow(presegmentation) > 0)
                for (presegRow in 1:nrow(presegmentation)) {
                  for (startToEnd in presegmentation[presegRow,2]:presegmentation[presegRow,2]) {
                    bedRow = bedFileForCluster[startToEnd,]
                    whichAreInside = which(bAlleleFreqsTumor[,1] == bedRow[,1] & as.numeric(bAlleleFreqsTumor[,2]) >= as.numeric(bedRow[,2]) & as.numeric(bAlleleFreqsTumor[,3]) <= as.numeric(bedRow[,3]))
                    notAllowedChromosomesAutosomesCNVs = c(notAllowedChromosomesAutosomesCNVs, whichAreInside)
                  }
                }
              
              allowedChromosomesAutosomesOnlySNV = setdiff(allowedChromosomesAutosomesOnlySNV, which(pvaluesSmall > 0.05 | pvaluesMed < 0.05))
              
              coordsIncludedAtFirst = unique(union(allowedChromosomesAutosomesOnlySNV, setdiff(1:nrow(bAlleleFreqsTumor), union(allowedChromosomesAutosomesOnlySNV, which(bAlleleFreqsTumor[,1] %in% ifelse(genderOfSamples[germline_sample_no] == "M",
                                                                                                                                                                                                            c("chrY"),
                                                                                                                                                                                                            c("chrX", "chrY")))))))
              if (length(coordsIncludedAtFirst) > 0) {
                bAlleleFreqsTumorToy = bAlleleFreqsTumor[coordsIncludedAtFirst,,drop=F]
                closestBedRegionsToy = closestBedRegions[coordsIncludedAtFirst]
                overdispersionNormalToy = overdispersionNormal[coordsIncludedAtFirst]
                overdispersionTumorToy = overdispersionTumor[coordsIncludedAtFirst]
              } else {
                noBAF = T
              }
            } else {
              bAlleleFreqsTumorToy = bAlleleFreqsTumor
              closestBedRegionsToy = closestBedRegions
              overdispersionNormalToy = overdispersionNormal
              overdispersionTumorToy = overdispersionTumor
            }
            if (!noBAF) {
              whichPUsedForAllStates = whichPTouseNew(local_purities, local_majorBAF, local_minorBAF, multiplierOfSNVsDueToMapping, degreeOfRoughness)
              lengthOfPositions <- sapply(1:length(whichPUsedForAllStates), function(i) {length(whichPUsedForAllStates[[i]])})
              expectedBAFsUnique = sort(unique(unlist(whichPUsedForAllStates)))
              whereThisBAFIsUsed <- sapply(1:length(expectedBAFsUnique),function(x) NULL)
              for (BAF in 1:length(expectedBAFsUnique)) {
                for (state in 1:length(whichPUsedForAllStates)) {
                  if (expectedBAFsUnique[BAF] %in% (whichPUsedForAllStates[[state]])) {
                    whereThisBAFIsUsed[[BAF]] = c(whereThisBAFIsUsed[[BAF]], state)
                  }
                }
              }
              matrixOfBAFLikeliks = foreach (i = 1:nrow(bAlleleFreqsTumorToy), .combine='rbind') %dopar% {
                altAlleleDepth <- as.numeric(bAlleleFreqsTumorToy[i,5])
                overallDepth <- round(as.numeric(bAlleleFreqsTumorToy[i,6]))
                altAlleleDepth = round(altAlleleDepth * overallDepth)
                overdispersionValue = overdispersionTumorToy[i]
                
                if (closestBedRegionsToy[i] != 0) {
                  numberOfAssignedPositions = numberOfAssignedPositions + 1
                  pList = fillInPList(altAlleleDepth, overallDepth,  expectedBAFsUnique, overdispersionValue)
                  vecOfLikeliks <- rep(0, length(local_cn_states))
                  for (j in 1:length(expectedBAFsUnique)) {
                    statesBAFIsUsed = whereThisBAFIsUsed[[j]]
                    likelihood = pList[[as.character(expectedBAFsUnique[j])]]
                    vecOfLikeliks[statesBAFIsUsed] = vecOfLikeliks[statesBAFIsUsed] + likelihood
                  }
                  -2 * log(vecOfLikeliks / lengthOfPositions)
                } else {
                  return(rep(0, length(local_cn_states)))
                }
              }
              #for (i in 1:nrow(bAlleleFreqsTumorToy)) {
              #  closestBedRegion = closestBedRegionsToy[i]
              #  if (!is.na(closestBedRegion))
              #    if (closestBedRegion != 0)
              #      matrix_of_likeliks[closestBedRegion,] = matrix_of_likeliks[closestBedRegion,] + matrixOfBAFLikeliks[i,]
              #}
            }
          }
          print("Finished BAF calculation")
          print(Sys.time())
        }
        # Reduce probability of unrealistic state = only super strong evidence is required to go for them
        
        
        
        sizesOfPointsFromLocalSds <- 0.5 / localSds 
        if (sampleInOfftarget) {
          matrixOfLogFoldOffCorrectedExtraSmallValues <- matrixOfLogFoldOff
          matrixOfLogFoldOffCorrectedExtraSmallValues[which(matrixOfLogFoldOffCorrectedExtraSmallValues < log2(min(local_cn_states) / 2))] = log2(min(local_cn_states) / 2)
          matrixOfLogFoldOffCorrectedExtraSmallValues[which(matrixOfLogFoldOffCorrectedExtraSmallValues > log2(max(local_cn_states) / 2))] = log2(max(local_cn_states) / 2)
          #matrix_of_likeliks_off <- form_matrix_of_likeliks_one_sample(1, matrixOfLogFoldOffCorrectedExtraSmallValues[,sam_no_off], localSdsOff, log2(local_cn_states/2))
          #if (genderOfSamples[germline_sample_no] == "M") {
          #  matrix_of_likeliks_off[which(bedFileForClusterOff[,1] %in% c("chrX","chrY")),] = form_matrix_of_likeliks_one_sample(
          #    1, matrixOfLogFoldOffCorrectedExtraSmallValues[which(bedFileForClusterOff[,1] %in% c("chrX","chrY")),sam_no_off], 
          #    localSdsOff[which(bedFileForClusterOff[,1] %in% c("chrX","chrY"))], log2((1 - local_purities) + local_purities * local_copy_numbers_used_major))
          #}
          
          
          
          #globalMatrOfLikeliks <- rbind(matrix_of_likeliks, matrix_of_likeliks_off)
          globalBed <- rbind(bedFileForCluster, bedFileForClusterOff)
          sizesOfPointsFromLocalSdsOff <- 0.5 / localSdsOff
          vecOfOrder = order(globalBed[,1], as.numeric(globalBed[,2]))
          globalSizesOfPoints <- c(sizesOfPointsFromLocalSds, sizesOfPointsFromLocalSdsOff)[vecOfOrder]
          #globalMatrOfLikeliks <- globalMatrOfLikeliks[vecOfOrder,]
          globalLogFold <- c( matrixOfLogFold[,sam_no], matrixOfLogFoldOff[,sam_no_off])[vecOfOrder]
          globalSds <-  c(localSds, localSdsOff)[vecOfOrder]
          globalBed <- globalBed[vecOfOrder,]
        }
        print(paste("Block before CNV detection finished", Sys.time()))
        
        #fileNameWithGermlineVars <- paste0(opt$out, "/normal/", germline_sample_name, "/", germline_sample_name, "_cnvs.tsv")
        #if (sampleInOfftarget) {
        #  coordsToMakeNull = returnCoordsThatNeedToBeNull(globalBed, fileNameWithGermlineVars)
        #  if (length(coordsToMakeNull) > 0)
        #    globalMatrOfLikeliks[coordsToMakeNull,] = 0
        #} else {
        #  coordsToMakeNull = returnCoordsThatNeedToBeNull(bedFileForCluster, fileNameWithGermlineVars)
        #  if (length(coordsToMakeNull) > 0)
        #    matrix_of_likeliks[coordsToMakeNull,] = 0
        #}
        
        
        
        found_CNVs_total <- matrix(0, nrow=0, ncol=14)
        colnames(found_CNVs_total) <- c("#chr", "start", "end", "major_CN_allele", "minor_CN_allele", "tumor_clonality", "CN_change", "loglikelihood", "median_loglikelihood", "number_of_regions", 
                                        "major_CN_allele2",  "minor_CN_allele2", "tumor_clonality2", "genes")
        allDetectedPurities = c()
        
        
        
        for (l in 1:length(left_borders)) {
          
          chrom = names(left_borders)[l]
          print(paste(chrom, Sys.time()))
          
          
          whichAreFromChr = which(bedFileForCluster[,1] == chrom)
          if (genderOfSamples[germline_sample_no]== "M" & chrom %in% c("chrX", "chrY")) {
            chrMatrixOfLikeliksOn = form_matrix_of_likeliks_one_sample(1,  matrixOfLogFoldCorrectedSmall[whichAreFromChr,sam_no], localSds[whichAreFromChr], log2((1 - local_purities) + local_purities * local_copy_numbers_used_major))
            if (sampleInOfftarget) chrMatrixOfLikeliksOff = form_matrix_of_likeliks_one_sample(1, matrixOfLogFoldOffCorrectedExtraSmallValues[which(bedFileForClusterOff[,1] == chrom),sam_no_off], localSdsOff[which(bedFileForClusterOff[,1] == chrom)], log2((1 - local_purities) + local_purities * local_copy_numbers_used_major))
          } else {
            chrMatrixOfLikeliksOn = form_matrix_of_likeliks_one_sample(1,  matrixOfLogFoldCorrectedSmall[whichAreFromChr,sam_no], localSds[whichAreFromChr], log2(local_cn_states/2))
            if (sampleInOfftarget) chrMatrixOfLikeliksOff = form_matrix_of_likeliks_one_sample(1, matrixOfLogFoldOffCorrectedExtraSmallValues[which(bedFileForClusterOff[,1] == chrom),sam_no_off], localSdsOff[which(bedFileForClusterOff[,1] == chrom)], log2(local_cn_states/2))
          }
          if (!is.null(nrow(chrMatrixOfLikeliksOn)))
            if (!is.null(matrixOfBAFLikeliks) & nrow(chrMatrixOfLikeliksOn) > 0) {
              closestBedRegionsToyChr = closestBedRegionsToy[which(bAlleleFreqsTumorToy[,1] == chrom)]
              matrixOfBAFsFromChr = matrixOfBAFLikeliks[which(bAlleleFreqsTumorToy[,1] == chrom),,drop=F]
              closestBedRegionsToyChrMapped = sapply(1:length(closestBedRegionsToyChr), function(i) {which(whichAreFromChr == closestBedRegionsToyChr[i])})
              for (mappedPosBaf in 1:length(closestBedRegionsToyChrMapped)) {
                if (length(closestBedRegionsToyChrMapped[[mappedPosBaf]]) > 0) {
                  chrMatrixOfLikeliksOn[closestBedRegionsToyChrMapped[[mappedPosBaf]],] = chrMatrixOfLikeliksOn[closestBedRegionsToyChrMapped[[mappedPosBaf]],] + matrixOfBAFsFromChr[mappedPosBaf,]
                }
              }
            }
          
          
          resultMatrixOfLikelihoods <- matrix(0, ncol=length(local_cn_states), nrow=0)
          resultBedForOrdering <- matrix(0, ncol=ncol(bedFileForCluster), nrow=0)
          if (!is.null(nrow(chrMatrixOfLikeliksOn)))
            if (nrow(chrMatrixOfLikeliksOn) > 0) {
              resultMatrixOfLikelihoods = rbind(resultMatrixOfLikelihoods, chrMatrixOfLikeliksOn)
              resultBedForOrdering = rbind(resultBedForOrdering, bedFileForCluster[which(bedFileForCluster[,1] == chrom),])
            }
          if (sampleInOfftarget) {
            if (nrow(chrMatrixOfLikeliksOff) > 0) {
              resultMatrixOfLikelihoods = rbind(resultMatrixOfLikelihoods, chrMatrixOfLikeliksOff)
              resultBedForOrdering = rbind(resultBedForOrdering, bedFileForClusterOff[which(bedFileForClusterOff[,1] == chrom),])
            }
          }
          resultMatrixOfLikelihoods = resultMatrixOfLikelihoods[order(resultBedForOrdering[,1], resultBedForOrdering[,2]),]
          resultBedForOrdering = resultBedForOrdering[order(resultBedForOrdering[,1], resultBedForOrdering[,2]),]
          
          
          
          if (chrom == "chrY" & genderOfSamples[germline_sample_no] == "F") {
            next
          } else if (genderOfSamples[germline_sample_no] == "M") {
            if (chrom %in% c("chrY","chrX")) {
              initial_state = 2
            } else {
              initial_state = 1
            }
          } else {
            initial_state = 1
          }
          
          
          
          start = left_borders[[l]]
          end = right_borders[[l]]
          for (k in 1:2) {
            
            output_of_plots <-  paste0(folder_name, sample_name)
            which_to_allow <- "NA"
            if (sampleInOfftarget) {
              if (k == 1) {
                which_to_allow = which(globalBed[,1] == chrom & as.numeric(globalBed[,2]) <= as.numeric(start) )
                which_to_allow_chr = which(as.numeric(resultBedForOrdering[,2]) <= as.numeric(start) )
              } else {
                which_to_allow = which(globalBed[,1] == chrom & as.numeric(globalBed[,2]) >= as.numeric(end) )
                which_to_allow_chr = which(as.numeric(resultBedForOrdering[,2]) >= as.numeric(end) )
              }
              toyBedFile = globalBed[which_to_allow,]
              
              toyMatrixOfLikeliks = resultMatrixOfLikelihoods[which_to_allow_chr,]
              toyLogFoldChange = globalLogFold[which_to_allow]
              
              toySds <- globalSds[which_to_allow]
              
              toySizesOfPointsFromLocalSds = c(sizesOfPointsFromLocalSds, sizesOfPointsFromLocalSdsOff)[vecOfOrder]
              toySizesOfPointsFromLocalSds = toySizesOfPointsFromLocalSds[which_to_allow]
            } else {
              if (k == 1) {
                which_to_allow = which(bedFileForCluster[,1] == chrom & bedFileForCluster[,2] <= start )
                which_to_allow_chr = which(as.numeric(resultBedForOrdering[,2]) <= as.numeric(start) )
              } else {
                which_to_allow = which(bedFileForCluster[,1] == chrom & bedFileForCluster[,2] >= end )
                which_to_allow_chr = which(as.numeric(resultBedForOrdering[,2]) >= as.numeric(end) )
              }
              toyMatrixOfLikeliks = resultMatrixOfLikelihoods[which_to_allow_chr,]
              toyBedFile = bedFileForCluster[which_to_allow,]
              toyLogFoldChange = matrixOfLogFold[which_to_allow, sam_no]
              
              
              toySds <- localSds[which_to_allow]
              
              toySizesOfPointsFromLocalSds = sizesOfPointsFromLocalSds[which_to_allow]
            }
            blocked_states = c()
            if (length(toyLogFoldChange) > opt$lengthS) {
              arrayOfMediansOfToyLogFold = runmed(toyLogFoldChange, round(opt$lengthS))
              if (!finalIteration) {
                diffsFromCoverage <- sapply(1:length(local_cn_states), function(i) {min(abs(log2(local_cn_states[i] / local_cn_states[initial_state]) - (arrayOfMediansOfToyLogFold)))})
                blocked_states = c(setdiff(c(1,2), initial_state),
                                   which(diffsFromCoverage > log2(1.05)))
              } else {
                blocked_states = c(setdiff(c(1,2), initial_state),
                                   which(log2(local_cn_states / local_cn_states[initial_state]) < min(arrayOfMediansOfToyLogFold) - 0.1 | log2(local_cn_states / local_cn_states[initial_state]) > max(arrayOfMediansOfToyLogFold) + 0.1))
              }
              if (genderOfSamples[germline_sample_no] == "M" & chrom %in% c("chrX","chrY")) {
                blocked_states = c(blocked_states, which(local_minorBAF != 0))
              }
              blocked_states = unique(blocked_states)
              if (initial_state %in% blocked_states) {
                blocked_states = blocked_states[-which(blocked_states == initial_state)]
              }
              if (length(local_cn_states) - length(blocked_states) == 1) {
                blocked_states <- c()
              }
            }
            if (!finalIteration) {
              datasetOfPuritiesCopiesSimplified = cbind(local_cn_states, local_majorBAF / (local_minorBAF + local_majorBAF))
              losingPrecision = which(duplicated(round(20 * datasetOfPuritiesCopiesSimplified, digits=0) / 20))
              blocked_states = unique(union(blocked_states, losingPrecision))
            }
            
            
            # BLOCK WITH PENALTIES
            copy_numbers_for_penalties = 3 - (local_copy_numbers_used_major + local_copy_numbers_used_minor)
            copy_numbers_for_penalties[which(copy_numbers_for_penalties > 0)] = 0
            penalties = penaltyForHigherCN * abs(copy_numbers_for_penalties)
            toyMatrixOfLikeliks = sweep(toyMatrixOfLikeliks, 2, abs(copy_numbers_for_penalties) * penaltyForHigherCNoneTile, FUN="+")
            #if (length(blocked_states) > 0) toyMatrixOfLikeliks[,blocked_states] = toyMatrixOfLikeliks[,blocked_states] + 10
            whichAreUnrealistic <- which((local_majorBAF == 0 & local_purities < 0.6) | local_purities_second > 10**-4)
            penalties[whichAreUnrealistic] = penalties[whichAreUnrealistic] + 20
            if (length(local_cn_states) - length(blocked_states) > 1) {
              found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, 1, blocked_states, penalties))
            } else {
              found_CNVs = matrix(0, nrow=0, ncol=10)
            }
            
            # BAFs from this chromosome
            if ((opt$filterStep == 1 & !finalIteration) | opt$filterStep >= 2) {
              if (frameworkDataTypes == "covdepthBAF" & !is.null(overdispersionNormal) & nrow(found_CNVs) > 0 & germline_sample_name %in% normalNames) {
                bafsFromThisChr = which(bAlleleFreqsNormal[,1] == chrom)
                listOfCNVsThatDoNotPass = returnListOfCNVsThatDoNotPass(found_CNVs, bafDeviationsForComparison, multiplierOfSNVsDueToMapping, bAlleleFreqsNormal[bafsFromThisChr,], bAlleleFreqsTumor[bafsFromThisChr,], 
                                                                        clonalityForChecking, local_purities, local_cn_states, toyBedFile,
                                                                        overdispersionNormal[bafsFromThisChr],
                                                                        overdispersionTumor[bafsFromThisChr],
                                                                        pvalueShift,
                                                                        toyLogFoldChange,
                                                                        median(sdsOfProbes) * (sdsOfSomaticSamples[sam_no]),
                                                                        ifelse(sampleInOfftarget, median(sdsOfProbesOff) * (sdsOfSomaticSamplesOff[sam_no_off]), -1)
                ) 
                if (length(listOfCNVsThatDoNotPass) > 0)
                  found_CNVs = found_CNVs[-listOfCNVsThatDoNotPass,,drop=F]
              }
            }
            
            
            if (nrow(found_CNVs) > 0){# & !chrom %in% c("chrX", "chrY", "X", "Y")) {
              likeliksFoundCNVsVsPurities <- matrix(0,nrow=nrow(found_CNVs), ncol=length(uniqueLocalPurities) * length(uniqueLocalPuritiesSecond))
              for (m in 1:length(uniqueLocalPurities)) {
                localPurityCurrent = uniqueLocalPurities[m]
                for (snd in 1:length(uniqueLocalPuritiesSecond)) {
                  localPurityCurrentSnd = uniqueLocalPuritiesSecond[snd]
                  for (q in 1:nrow(found_CNVs)) {
                    startOfCNV <- found_CNVs[q,2]
                    endOfCNV <- found_CNVs[q,3]
                    if (endOfCNV - startOfCNV > 3) { 
                      
                      likeliksFoundCNVsVsPurities[q, (snd - 1) * m + m] = min(apply(toyMatrixOfLikeliks[(startOfCNV + 1):(endOfCNV - 1),which(local_purities == localPurityCurrent & local_purities_second == localPurityCurrentSnd),drop=F], 2, sum)
                                                                              + penaltyForHigherCN * abs(copy_numbers_for_penalties[which(local_purities == localPurityCurrent & local_purities_second == localPurityCurrentSnd)]))
                    }
                  }
                }
              }
              likeliksFoundCNVsVsPuritiesGlobal = rbind(likeliksFoundCNVsVsPuritiesGlobal, likeliksFoundCNVsVsPurities)
            }
            if (opt$visulizationIGV) {
              ### IGV PLOTTING
              if(opt$debug) {
                print(paste("Start of IGV plotting", Sys.time()))
              }
              if (finalIteration) {
                outputFileNameCNVs <- paste0(folder_name, sample_name, "/", sample_name, "_cnvs.seg")
                outputFileNameDots <- paste0(folder_name, sample_name, "/", sample_name, "_cov.seg")
                if (genderOfSamples[germline_sample_no] == "M") {
                  reverseFunctionUsedToTransform = function(x, chroms) {return((2 ** (x + ifelse(chrom %in% c("chrX", "chrY"), 0, 1))))}
                } else {
                  reverseFunctionUsedToTransform = function(x, chroms) {return(2 ** (x + 1))}
                }
                outputSegmentsAndDotsFromListOfCNVs(toyBedFile, found_CNVs, start, end, outputFileNameCNVs, 
                                                    outputFileNameDots, sample_name, toyLogFoldChange, reverseFunctionUsedToTransform, local_cn_states, toySds)
              }
              if(opt$debug) {
                print(paste("End of IGV plotting", Sys.time()))
              }
            }
            if (length(local_cn_states) - length(blocked_states) <= 1) {
              next
            }       
            ### END OF IGV PLOTTING
            
            
            if (nrow(found_CNVs) == 0 & length(which_to_allow) > 1) {
              found_CNVs = matrix(c(-1000, 1, length(which_to_allow), 1), nrow=1)
            }
            
            
            
            if (nrow(found_CNVs) > 0) {
              medianLikelihoods <- round(-1 * sapply(1:nrow(found_CNVs), function(i) {median(toyMatrixOfLikeliks[found_CNVs[i,2]:found_CNVs[i,3], found_CNVs[i,4]] - 
                                                                                               toyMatrixOfLikeliks[found_CNVs[i,2]:found_CNVs[i,3], initial_state])}), 3)
              cnvsToWriteOut <- plotFoundCNVsNew(sam_no, found_CNVs, toyLogFoldChange, toyBedFile, output_of_plots, chrom, 
                                                 local_cn_states, local_copy_numbers_used_major, local_copy_numbers_used_minor, local_purities, 
                                                 local_copy_numbers_used_major_second, local_copy_numbers_used_minor_second, local_purities_second,
                                                 toySizesOfPointsFromLocalSds, plottingOfPNGs, medianLikelihoods)
              if (found_CNVs[1,1] != -1000) {
                found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
              }
            }
          }
        }
        ### Sometimes false positive CNV calls can be caused by deviation in coverage in normal - we correct it
        if (!finalIteration) {
          if (nrow(found_CNVs_total) > 0) {
            if (sampleInOfftarget) {
              shiftsForCoverageInsideCNVs <- findDeviationInNormalCoverage(germline_sample_name, tumor_sample_name, found_CNVs_total, bedFileForCluster, tmpNormal,
                                                                           bedFileForClusterOff, tmpNormalOff)
            } else { 
              shiftsForCoverageInsideCNVs <- findDeviationInNormalCoverage(germline_sample_name, tumor_sample_name, found_CNVs_total, bedFileForCluster, tmpNormal)
            }
            
            for (m in 1:nrow(shiftsForCoverageInsideCNVs)) {
              cnvSpecificShift = log2(shiftsForCoverageInsideCNVs[m,1])
              chrom = found_CNVs_total[m,1]
              start = as.numeric(found_CNVs_total[m,2])
              end = as.numeric(found_CNVs_total[m,3])
              coordsInBedFileOn = which(bedFileForCluster[,1] == chrom & as.numeric(bedFileForCluster[,2]) >= start & as.numeric(bedFileForCluster[,3]) <= end)
              matrixOfLogFold[coordsInBedFileOn,sam_no] = matrixOfLogFold[coordsInBedFileOn,sam_no] + cnvSpecificShift
              if (sampleInOfftarget) {
                coordsInBedFileOff = which(bedFileForClusterOff[,1] == chrom & as.numeric(bedFileForClusterOff[,2]) >= start & as.numeric(bedFileForClusterOff[,3]) <= end)
                matrixOfLogFoldOff[coordsInBedFileOff,sam_no_off] = matrixOfLogFoldOff[coordsInBedFileOff,sam_no_off] + cnvSpecificShift
              }
            }
            if (length(which(shiftsForCoverageInsideCNVs[,2] == 0)) > 0) {
              found_CNVs_total = found_CNVs_total[-which(shiftsForCoverageInsideCNVs[,2] == 0),,drop=F]
              likeliksFoundCNVsVsPuritiesGlobal = likeliksFoundCNVsVsPuritiesGlobal[-which(shiftsForCoverageInsideCNVs[,2] == 0),,drop=F]# & !found_CNVs_total[,1] %in% c("chrX","chrY")),,drop=F]
            }
          }
        }
        
        if (finalIteration == T) {
          if (nrow(found_CNVs_total) > 0){
            # HOMOZYGOUSITY FILTER - IF THERE ARE TOO MANY HOMOZYGOUS, WE REMOVE THEM
            if (length( which(as.numeric(found_CNVs_total[,4]) == 0 & as.numeric(found_CNVs_total[,9]) < 10 & as.numeric(found_CNVs_total[,6]) > 0.5) ) > 5) {
              print("Short (<10 regions) homozygous deletions will be filtered out due to high percentage of technical artifacts in such CNVs")
              print(found_CNVs_total[which((as.numeric(found_CNVs_total[,4]) == 0 & as.numeric(found_CNVs_total[,9]) < 10 & as.numeric(found_CNVs_total[,6]) > 0.5)),,drop=F])
              found_CNVs_total = found_CNVs_total[which(!(as.numeric(found_CNVs_total[,4]) == 0 & as.numeric(found_CNVs_total[,9]) < 10 & as.numeric(found_CNVs_total[,6]) > 0.5)),,drop=F]
            }
            makeBarplot(allPotentialPurities, found_CNVs_total, sample_name)
            plotChromosomalLevelInstabs(found_CNVs_total, left_borders, right_borders, ends_of_chroms, genderOfSamples[germline_sample_no], sample_name)
          }
          
          break}
        
        
        
        
        
        
        if (nrow(likeliksFoundCNVsVsPuritiesGlobal) > 0 & !finalIteration) {
          # Again - matrix of clonality
          matrixOfClonality = matrix(0, nrow=length(uniqueLocalPurities), ncol=length(uniqueLocalPurities))
          colnames(matrixOfClonality) = uniqueLocalPurities
          rownames(matrixOfClonality) = uniqueLocalPurities
          for (m in 1:length(uniqueLocalPurities)) {
            for (q in 1:length(uniqueLocalPurities)) {
              for (r in 1:nrow(likeliksFoundCNVsVsPuritiesGlobal)) {
                matrixOfClonality[m,q] = matrixOfClonality[m,q] + min(likeliksFoundCNVsVsPuritiesGlobal[r,m], likeliksFoundCNVsVsPuritiesGlobal[r,q])
              }
            }
          }
          minPointToNormalise = quantile(matrixOfClonality, 0.25)
          matrixOfClonalityForPlotting = matrixOfClonality
          matrixOfClonalityForPlotting[which(matrixOfClonalityForPlotting > minPointToNormalise)] = minPointToNormalise
          matrixOfClonalityForPlotting[upper.tri(matrixOfClonalityForPlotting)] <- NA
          hmcols<-colorRampPalette(c("blue","white","red"))(256)
          if (!finalIteration) {
            do.call(file.remove, list(list.files(paste0(folder_name, sample_name), full.names = TRUE)))
          }
          png(filename = paste0(sample_name, "_clonality.png"),
              width = 640, height = 640)
          heatmap((matrixOfClonalityForPlotting), scale="none", Rowv = NA, Colv = NA, col=hmcols, main=sample_name)
          dev.off()
          
          
          # true clonal structure
          maxNumOfClones <- min(5, length(uniqueLocalPurities))
          localPurityStates = 1:length(uniqueLocalPurities)
          #hundredPercentPurity = which(uniqueLocalPurities == 1)
          resultBestCombination = 0
          minResult = 10**100
          for (m in 2:maxNumOfClones) {
            print(paste("Investigating combinations of", m, "clones"))
            combinationsOfPurities <- combn(localPurityStates, m)
            #indicesOfPuritiesWithMax = apply(combinationsOfPurities, 2, function(x) {sum(which(x == hundredPercentPurity))})
            #combinationsOfPurities = combinationsOfPurities[,which(indicesOfPuritiesWithMax > 0),drop=F]
            
            bestCombination = 1
            
            for (q in 1:ncol(combinationsOfPurities)) {
              minResultForCombination = as.numeric(opt$clonePenalty) * m
              positionsWithSpecifiedPurities = combinationsOfPurities[,q]
              if (length(intersect(
                uniqueLocalPurities[combinationsOfPurities[,q]], uniqueLocalPuritiesSecond)) > 0) {
                whichToIncludeFromSecond <- which(uniqueLocalPuritiesSecond %in% uniqueLocalPurities[combinationsOfPurities[,q]])
                for (sndPurity in whichToIncludeFromSecond) {
                  positionsWithSpecifiedPurities = c(positionsWithSpecifiedPurities, (sndPurity - 1) * combinationsOfPurities[,q] + combinationsOfPurities[,q])
                }
                positionsWithSpecifiedPurities = unique(positionsWithSpecifiedPurities)
              }
              for (r in 1:nrow(likeliksFoundCNVsVsPuritiesGlobal)) {
                minResultForCombination = minResultForCombination + min(
                  likeliksFoundCNVsVsPuritiesGlobal[r,positionsWithSpecifiedPurities])
              }
              if (minResult > minResultForCombination) {
                bestCombination = q
                resultBestCombination = combinationsOfPurities[,bestCombination]
                minResult = minResultForCombination
              }
            }
          }
          print("Clonal sturcture inferred.")
          clonalBestPurities = uniqueLocalPurities[resultBestCombination]
          if (length(clonalBestPurities) == 0) {
            clonalBestPurities = c(0, 1)
          }
          print(clonalBestPurities)
          
        } else {
          print("No high quality CNVs found in this sample for finding clonality.")
          clonalBestPurities = c(0, 1)
        }
        finalIteration = T
      }
      addressOfPlot = paste0(sample_name, "_CNAs_plot.png")
      if (finalIteration & frameworkDataTypes == "covdepthBAF") {
        if (sampleInOfftarget) {
          plotLikelihoodLandscape(datasetOfPuritiesCopiesForFinalIteration, addressOfPlot, found_CNVs_total, globalBed, matrixOfBAFLikeliks, bAlleleFreqsTumor, 
                                  coordsIncludedAtFirst, globalLogFold, local_purities, local_majorBAF, local_minorBAF, left_borders, right_borders, ends_of_chroms,
                                  local_copy_numbers_used_major_second, local_copy_numbers_used_minor_second, local_purities_second)
        } else {
          plotLikelihoodLandscape(datasetOfPuritiesCopiesForFinalIteration, addressOfPlot, found_CNVs_total, bedFileForCluster, matrixOfBAFLikeliks, bAlleleFreqsTumor, 
                                  coordsIncludedAtFirst,matrixOfLogFold[,sam_no],local_purities, local_majorBAF, local_minorBAF, left_borders, right_borders, ends_of_chroms,
                                  local_copy_numbers_used_major_second, local_copy_numbers_used_minor_second, local_purities_second)
        }
      }
      
      ### STAT TESTS
      CIsOnTarget = matrix(NA, nrow=nrow(found_CNVs_total), ncol=3)
      CIsOnTargetOff = matrix(NA, nrow=nrow(found_CNVs_total), ncol=3)
      BAFsignature = matrix(NA, nrow=nrow(found_CNVs_total), ncol=3)
      snvUpperAndBottom = matrix(NA, nrow=nrow(found_CNVs_total), ncol=2)
      overallPvalues = matrix(NA, nrow=nrow(found_CNVs_total), ncol=1)
      QC_value = "NA"
      if (nrow(found_CNVs_total) > 0){
        for (i in 1:nrow(found_CNVs_total)) {
          defaultCN = 2
          if (found_CNVs_total[i,1] %in% c("chrX", "chrY") & genderOfSamples[germline_sample_no] == "M") {
            defaultCN = 1
          }
          pvalsSeparateTests = c(NA, NA, NA)
          onTargetCoords <- which(bedFileForCluster[,1] == found_CNVs_total[i,1] & as.numeric(bedFileForCluster[,2]) >= as.numeric(found_CNVs_total[i,2]) & as.numeric(bedFileForCluster[,3]) <= as.numeric(found_CNVs_total[i,3]))
          if (length(onTargetCoords) > 1) {
            tumorValue <- median(log2(tmpTumor[onTargetCoords, tumor_sample_no]))
            if (found_CNVs_total[i,1] %in% c("chrX", "chrY")) {
              samplesToUse = which(genderOfSamples == genderOfSamples[germline_sample_no])
            } else {
              samplesToUse = 1:ncol(tmpNormal)
            }
            if (length(samplesToUse) > 2) {
              normalValues <- apply((log2(tmpNormal[onTargetCoords,samplesToUse, drop=F])), 2, median)
              sdOfNormals <- sd(normalValues) * sqrt(matrixWithSds[1, sam_no] / median(matrixWithSds[2, ]))
              currentCI = c((tumorValue), (tumorValue + qnorm(0.99) * sdOfNormals), (tumorValue + qnorm(0.01) * sdOfNormals))
              currentCI = 2 ** (currentCI - median(log2(tmpNormal[onTargetCoords, germline_sample_no])))
              CIsOnTarget[i,] = round(currentCI * defaultCN, 2)
            }
            pvalsSeparateTests[1] = 2 * pt(-abs(   
              tumorValue - median(log2(tmpNormal[onTargetCoords, germline_sample_no]))
            ) / sdOfNormals, df=length(samplesToUse)) 
          }
          if (sampleInOfftarget) {
            offTargetCoords <- which(bedFileForClusterOff[,1] == found_CNVs_total[i,1] & as.numeric(bedFileForClusterOff[,2]) >= as.numeric(found_CNVs_total[i,2]) & as.numeric(bedFileForClusterOff[,3]) <= as.numeric(found_CNVs_total[i,3]))
            if (length(offTargetCoords) > 1) {
              tumorValueOff <- median(log2(tmpTumorOff[offTargetCoords, tumor_sample_no_off]))
              if (found_CNVs_total[i,1] %in% c("chrX", "chrY")) {
                samplesToUseOn = which(genderOfSamples == genderOfSamples[germline_sample_no])
                samplesToUseOff = which(colnames(tmpNormalOff) %in% names(samplesToUseOn))
              } else {
                samplesToUseOff = 1:ncol(tmpNormalOff)
              }
              if (length(samplesToUseOff) > 2) {
                normalValuesOff <- apply(log2(tmpNormalOff[offTargetCoords,samplesToUseOff, drop=F]), 2, median)
                sdOfNormalsOff <- sd(normalValuesOff) * sqrt(matrixWithSdsOff[1, sam_no_off] / median(matrixWithSdsOff[2, ]))
                currentCIOff = c((tumorValueOff), (tumorValueOff + qnorm(0.99) * sdOfNormalsOff), (tumorValueOff + qnorm(0.01) * sdOfNormalsOff))
                currentCIOff = 2 ** (currentCIOff - median(log2(tmpNormalOff[offTargetCoords, which(colnames(tmpNormalOff) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2])])))
                CIsOnTargetOff[i,] = round(currentCIOff * defaultCN, 2)
              }
              pvalsSeparateTests[2] = 2 * pt( -abs(
                tumorValueOff - median(log2(tmpNormalOff[offTargetCoords, which(colnames(tmpNormalOff) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2])])))
                / sdOfNormalsOff, 
                df=length(samplesToUseOff)) 
            }
          }
          if (!is.null(bAlleleFreqsTumor) & !is.null(bAlleleFreqsNormal)) {
            BAFcoords <- which(bAlleleFreqsNormal[,1] == found_CNVs_total[i,1] & as.numeric(bAlleleFreqsNormal[,2]) >= as.numeric(found_CNVs_total[i,2]) & as.numeric(bAlleleFreqsNormal[,2]) <= as.numeric(found_CNVs_total[i,3]))
            particularAlleleBalance = median(as.numeric(bAlleleFreqsNormal[,5]))
            
            if (length(BAFcoords) > 0) {
              frequencies = as.numeric(bAlleleFreqsNormal[BAFcoords,5])
              frequenciesTumor = as.numeric(bAlleleFreqsTumor[BAFcoords,5])
              snvUpperAndBottom[i,1] = round(median(frequenciesTumor[which(frequenciesTumor < particularAlleleBalance)]),3)
              snvUpperAndBottom[i,2] = round(median(frequenciesTumor[which(frequenciesTumor > particularAlleleBalance)]),3)
              depthNormal = sum(as.numeric(bAlleleFreqsNormal[BAFcoords,6]))
              normalReversedValues = frequencies
              centralLocation = median(normalReversedValues)
              diffs = abs(normalReversedValues - centralLocation)
              normalReversedValues = centralLocation + diffs
              #normalReversedValues[which(normalReversedValues > particularAlleleBalance)] = max(0, 2 * particularAlleleBalance - normalReversedValues[which(normalReversedValues > particularAlleleBalance)])
              medianAFNormal = median(normalReversedValues)
              normalDepthRefAllele = round(medianAFNormal * depthNormal)
              
              depthTumor = sum(as.numeric(bAlleleFreqsTumor[BAFcoords,6]))
              tumorReversedValues = as.numeric(bAlleleFreqsTumor[BAFcoords,5])
              diffs = abs(tumorReversedValues - centralLocation)
              tumorReversedValues = centralLocation + diffs
              #tumorReversedValues[which(tumorReversedValues > particularAlleleBalance)] = max(0, 2 * particularAlleleBalance - tumorReversedValues[which(tumorReversedValues > particularAlleleBalance)])
              medianAFTumor = min(1, median(tumorReversedValues))
              tumorDepthRefAllele = round(medianAFTumor * depthTumor)
              
              BAFsignature[i,1] = paste(normalDepthRefAllele, "/", depthNormal)
              BAFsignature[i,2] = paste(tumorDepthRefAllele, "/", depthTumor)
              
              BAFsignature[i,3] = (fisher.test(matrix(c(normalDepthRefAllele, depthNormal - normalDepthRefAllele, tumorDepthRefAllele, depthTumor - tumorDepthRefAllele), nrow=2))$p.value)
              pvalsSeparateTests[3] = BAFsignature[i,3]
            }
            
          }
          pvalsSeparateTests = as.numeric(na.omit((pvalsSeparateTests)))
          overallPvalues[i] = pchisq((sum(log(min(1, pvalsSeparateTests + 10**-10)))*-2), df=length(pvalsSeparateTests)*2, lower.tail=F)
        }
        overallPvalues = p.adjust(overallPvalues, method="fdr")
        if (length(which(found_CNVs_total[,4] != found_CNVs_total[,5])) > 0) {
          QC_value = 2 * length(which(as.numeric(BAFsignature[which(found_CNVs_total[,5] != found_CNVs_total[,4]),3]) > 0.5)) / length(which(found_CNVs_total[,4] != found_CNVs_total[,5]))
        } else {
          QC_value = "NA"
        }
        BAFsignature[,3] = p.adjust(as.numeric(BAFsignature[,3]), method="fdr")
        BAFsignature[,3] = format(round(as.numeric(BAFsignature[,3]), 4), scientific = F)
        colnamesForFutureMatrix <- colnames(found_CNVs_total)
        tumor_cn_change = as.numeric(found_CNVs_total[,4]) + as.numeric(found_CNVs_total[,5])
        tumor_cn_state = round(as.numeric(found_CNVs_total[,7]), 2)
        states_for_report = sapply(1:nrow(found_CNVs_total), function(i) {if (tumor_cn_state[i] == 2) return("LOH"); if (tumor_cn_state[i] > 2) return("AMP"); if (tumor_cn_state[i] < 2) return("DEL")})
        colnamesForFutureMatrix = c(colnamesForFutureMatrix[1:3], "tumor_CN_change", "state", colnamesForFutureMatrix[4:length(colnamesForFutureMatrix)])
        found_CNVs_total = cbind(found_CNVs_total[,1:3,drop=F], tumor_cn_change, states_for_report, found_CNVs_total[,4:ncol(found_CNVs_total),drop=F],
                                 matrix(CIsOnTarget[,3:2,drop=F], ncol=2), matrix(CIsOnTargetOff[,3:2,drop=F], ncol=2), snvUpperAndBottom, BAFsignature[,3,drop=F], format(round(overallPvalues,5), scientific = F))
        #colnames(found_CNVs_total) = c(colnamesForFutureMatrix, c("Ontarget_RD", "Ontarget_RD_CI_lower", "Ontarget_RD_CI_upper", "Offtarget_RD", "Offtarget_RD_CI_lower", "Offtarget_RD_CI_upper", "BAF_Normal", "BAF_tumor", "BAF_pval"))
        colnames(found_CNVs_total) = c(colnamesForFutureMatrix, c("Ontarget_RD_CI_lower", "Ontarget_RD_CI_upper", "Offtarget_RD_CI_lower", "Offtarget_RD_CI_upper", "Lowmed_tumor_BAF", "Highmed_tumor_BAF",  "BAF_qval_fdr", "Overall_qvalue"))
      }
      if (length(pvalsForQC > 1)) {
        finalPValue <- 0
      } else {
        finalPValue = 0
      }
      fileToOut <- paste0(folder_name, sample_name, paste0("/CNAs_", sample_name, ".txt"))
      fileConn<-file(fileToOut)
      ploidyEst = round(2 + (2 * 2 ** median(matrixOfLogFold[which(!bedFileForCluster[,1] %in% c("chrX", "chrY")),sam_no]) - 2) / ifelse(max(local_purities) == 0, 1, max(local_purities)), 4)
      estimatedFDR = ifelse(QC_value == "NA", "NA", round(QC_value, 5))
      writeLines(c(paste0("##ANALYSISTYPE=CLINCNV_TUMOR_NORMAL_PAIR"), paste0( "##", clincnvVersion), paste0("##Analysis finished on: ", Sys.time()),  paste0("##estimated fdr: ", estimatedFDR), paste0("##gender of sample: ", genderOfSamples[germline_sample_no]), paste0("##ploidy: ", ploidyEst), paste0("##clonality by BAF (if != 1): ", paste(round(unique(local_purities), digits=3), collapse=";"), collapse = " ")), fileConn)      close(fileConn)
      
      if(opt$debug) {
        print(found_CNVs_total)
      }
      found_CNVs_total = found_CNVs_total[order(found_CNVs_total[,1], as.numeric(found_CNVs_total[,2])), , drop=F]
      if (nrow(found_CNVs_total) > 0) found_CNVs_total[which(found_CNVs_total[,15] == "0"), 13:15] = ""
      write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)	
      # For some additional analysis we need to provide areas free of CNVs
      if (sampleInOfftarget) {
        areasFreeOfCNVs <- returnAreasFreeOfCNVsForAdditionalAnalysis(found_CNVs_total, genderOfSamples[germline_sample_no], bedFileForCluster, bedFileForClusterOff)
      } else {
        areasFreeOfCNVs <- returnAreasFreeOfCNVsForAdditionalAnalysis(found_CNVs_total, genderOfSamples[germline_sample_no], bedFileForCluster)
      }
      if (nrow(areasFreeOfCNVs) > 0) {
        fileToOut <- paste0(folder_name, sample_name, paste0("/CNneutral_", sample_name, ".txt"))
        write.table(areasFreeOfCNVs, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)
      }
    }
  }
}

