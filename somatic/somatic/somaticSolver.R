library(mclust)

### PROCESSING OF SOMATIC VARIANTS
setwd(opt$folderWithScript)
source("./somatic/helpersSomatic.R",local=TRUE)

listOfValue <- formilngLogFoldChange(pairs, normal, tumor)
matrixOfLogFold <- listOfValue[[1]]
dictFromColumnToTumor <- listOfValue[[2]]

bordersOfChroms <- getBordersOfChromosomes(bedFile)
#sdsOfSomaticSamples <- sapply(1:ncol(matrixOfLogFold), function(i) {determineSDsOfSomaticSample(matrixOfLogFold[,i], bedFile)})
matrixWithSds <- findSDsOfSamples(pairs, normal, tumor, bedFile, bordersOfChroms)
sdsOfSomaticSamples <- matrixWithSds[4,]

sdsOfProbes <- sapply(1:nrow(matrixOfLogFold), function(i) {determineSDsOfSomaticProbe(matrixOfLogFold[i,], i)})

listOfVarianceAndMultiplicator <- esimtateVarianceFromSampleNoise(sdsOfSomaticSamples, 100000)
esimtatedVarianceFromSampleNoise <- listOfVarianceAndMultiplicator[[2]]
multiplicator <- listOfVarianceAndMultiplicator[[1]]

if (frameworkOff == "offtarget") {
  listOfValue <- formilngLogFoldChange(pairs, normalOff, tumorOff)
  matrixOfLogFoldOff =  listOfValue[[1]]
  dictFromColumnToTumorOff = listOfValue[[2]]
  bordersOfChroms <- getBordersOfChromosomes(bedFileOfftarget)
  matrixWithSdsOff <- findSDsOfSamples(pairs, normalOff, tumorOff, bedFileOfftarget, bordersOfChroms)
  sdsOfSomaticSamplesOff <- matrixWithSdsOff[4,]

  sdsOfProbesOff <- sapply(1:nrow(matrixOfLogFoldOff), function(i) {determineSDsOfSomaticProbe(matrixOfLogFoldOff[i,], i)})
  
  listOfVarianceAndMultiplicatorOff <- esimtateVarianceFromSampleNoise(sdsOfSomaticSamplesOff, 100000)
  esimtatedVarianceFromSampleNoiseOff <- listOfVarianceAndMultiplicatorOff[[2]]
  multiplicatorOff <- listOfVarianceAndMultiplicatorOff[[1]]
  
}

### FORMING MATRIX OF LIKELIHOODS
vect_of_t_likeliks <- fast_dt_list(ncol(matrixOfLogFold) - 1)




# cn_states <- c()
# purity <- seq(from=10, to=101, by=5) / 100
# for (pur in purity) {
#   cn_state <- 0 + 1 * pur * seq(from=0, to=20, by=1)
#   cn_states <- c(cn_states, cn_state)
# }
# cn_states <- unique(cn_states[-which(cn_states > 1.5 & cn_states < 2.8)])
# cn_state[which(cn_state < 0.01)] = 0.01


cn_states <- c()
copy_numbers = 0:15
purity <- seq(from=5, to=101, by=5) / 100
purities <- c()
copy_numbers_used <- c()
statesUsed <- c()

### DESCRIPTION OF STATES
# CNV - copy number change, 1 allele changed
# CNVboth - duplication when both alleles changed
# LOH - Loss of Heterozygosity
# normal - nothing changed comparing to normal genome
# CNVcomplex - not single allelic CNV

for (pur in purity) {
  for (cn in copy_numbers) {
    # VALUES NOT CN NEUTRAL
    if (cn != 2 & (!(cn == 0 & pur <= 0.5))) {
      cn_state <- (1 - pur) * 2 + pur * cn
      cn_states <- c(cn_states, cn_state)
      purities <- c(purities, pur)
      copy_numbers_used <- c(copy_numbers_used, cn)
      statesUsed <- c(statesUsed, "CNV")
    }
    # CN neutral changes
    if (cn == 2) {
      if (pur >= 0.1) {
        cn_states <- c(cn_states, 2)
        copy_numbers_used <- c(copy_numbers_used, 2)
        purities <- c(purities, pur)
        statesUsed <- c(statesUsed, "LOH")
      }
    }
    if (cn == 4 | cn == 6 | cn == 8) {
      cn_state <- (1 - pur) * 2 + pur * cn
      cn_states <- c(cn_states, cn_state)
      copy_numbers_used <- c(copy_numbers_used, cn)
      purities <- c(purities, pur)
      statesUsed <- c(statesUsed, "CNVboth")
    }
    if (cn >= 5 & cn <= 8) {
      cn_state <- (1 - pur) * 2 + pur * cn
      cn_states <- c(cn_states, cn_state)
      copy_numbers_used <- c(copy_numbers_used, cn)
      purities <- c(purities, pur)
      statesUsed <- c(statesUsed, "CNVcomplex")
    }
    if (cn == 3 | cn == 4) {
      cn_state <- (1 - pur) * 2 + pur * cn
      cn_states <- c(cn_states, cn_state)
      copy_numbers_used <- c(copy_numbers_used, cn)
      purities <- c(purities, pur)
      statesUsed <- c(statesUsed, "LOHDup")
    }
  }
}




cn_states[which(cn_states < 0.0001)] = 0.0001
matrixOfLogFold[which(matrixOfLogFold < log2(min(cn_states) / 2))] = log2(min(cn_states) / 2)
matrixOfLogFold[which(matrixOfLogFold > log2(max(cn_states) / 2))] = log2(max(cn_states) / 2)
if (frameworkOff == "offtarget") {
  matrixOfLogFoldOff[which(matrixOfLogFoldOff < log2(min(cn_states) / 2))] = log2(min(cn_states) / 2)
  matrixOfLogFoldOff[which(matrixOfLogFoldOff > log2(max(cn_states) / 2))] = log2(max(cn_states) / 2)
}

final_order <- order(cn_states)
cn_states <- cn_states[final_order]
copy_numbers_used = copy_numbers_used[final_order]
purities = purities[final_order]
statesUsed = statesUsed[final_order]
### ADDING NORMAL STATE
purities <- c(0, 0, purities)
copy_numbers_used <- c(2, 1, copy_numbers_used)
cn_states <- c(2, 1, cn_states)
statesUsed <- c("normal", "normal", statesUsed)






cnsLessThanTwo <- which(cn_states < 2)
cnsBiggerThanTwo <- which(cn_states > 2 & cn_states < 4)
cnsHighCopies <-  which(cn_states >= 4)

colfunc <- colorRampPalette(c("darkblue", "lightblue"))
colorsForLess <- colfunc(length(cnsLessThanTwo))

colfunc <- colorRampPalette(c("yellow", "brown"))
colorsForBigger <- colfunc(length(cnsBiggerThanTwo))

colfunc <- colorRampPalette(c("red", "darkred"))
colorsForHigh <- colfunc(length(cnsHighCopies))

colours <- c("darkgreen", colorsForLess, colorsForBigger, colorsForHigh)



# CORRECTION OF CNS!!!
# sdsOfNormals <- apply(normal, 1, sd)
# medianBaselineSD <- median(sdsOfNormals)
# 
# sdsPois <- c()
# for (i in 1:1000) {
#   samplePois <- rpois(lambda=i, n=10000)
#   sdsPois <- c(sdsPois, sd(samplePois / i))
# }
# diffs <- abs(sdsPois - medianBaselineSD)
# lambdaFromSimulation <- which.min(diffs)
# 
# normalCoverage <- rpois(1000000, lambda=lambdaFromSimulation)
# tumorCoverage <- rpois(1000000, lambda=lambdaFromSimulation)
# sd_to_normalise = sd(log2(tumorCoverage / normalCoverage))
# multipliersDueToLog <- c(1)
# for (state in 2:length(cn_states)) {
#   tumorCoverage <- rpois(1000000, lambda=0.5 * lambdaFromSimulation * cn_states[state])
#   sd_to_normalise_tumor = sd(log2(tumorCoverage / normalCoverage))
#   multipliersDueToLog <- c(multipliersDueToLog, sd_to_normalise_tumor / sd_to_normalise)
# }
# multipliersDueToLog[which(is.nan(multipliersDueToLog))] <- max(multipliersDueToLog[which(!is.nan(multipliersDueToLog))])
# 
# ### Make SDs for states equal
# for (state in unique(cn_states)) {
#   whichStatesAre = which(cn_states == state)
#   if (state != 2) {
#     multipliersDueToLog[whichStatesAre] = median(multipliersDueToLog[whichStatesAre])
#   } else {
#     multipliersDueToLog[whichStatesAre] = 1
#   }
# }





folder_name <- paste0(opt$out, "/somatic/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

allPotentialPurities <- unique(purities)
penaltyForHigherCN = 20
for (sam_no in 1:ncol(matrixOfLogFold)) {
  sample_name <- colnames(matrixOfLogFold)[sam_no]
  germline_sample_no = which(colnames(normal) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2])
  if (germline_sample_no == "F") next
  print(genderOfSamples[germline_sample_no])
  
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
  
  ### Multiplier due to different CN - calculation with covariance
  multipliersDueToLog = c()
  for (cn in cn_states) {
  multipliersDueToLog = c(multipliersDueToLog, returnMultiplierDueToLog(2,cn,matrixWithSds[1,sam_no], matrixWithSds[2,sam_no], matrixWithSds[3,sam_no]))
  }
  multipliersDueToLog = sqrt(multipliersDueToLog) 
  multipliersDueToLog = multipliersDueToLog / multipliersDueToLog[1]
  
  if (frameworkOff == "offtarget") {
    if (sample_name %in% colnames(matrixOfLogFoldOff)) {
      sam_no_off = which(colnames(matrixOfLogFoldOff) == sample_name)
      multipliersDueToLogOff = c()
      for (cn in cn_states) {
        multipliersDueToLogOff = c(multipliersDueToLogOff, returnMultiplierDueToLog(2,cn,matrixWithSdsOff[1,sam_no_off], matrixWithSdsOff[2,sam_no_off], matrixWithSdsOff[3,sam_no_off]))
      }
      multipliersDueToLogOff = sqrt(multipliersDueToLogOff) 
      multipliersDueToLogOff = multipliersDueToLogOff / multipliersDueToLogOff[1]
    }
  }
  
  
  matrixOfClonality = matrix(0, nrow=1, ncol=1)
  if (!dir.exists(paste0(folder_name, sample_name)) | (opt$reanalyseCohort == T)) {
    
    dir.create(paste0(folder_name, sample_name))
    setwd(paste0(folder_name, sample_name))
    
    copyNumbersInsideExpectedPurities = F
    
    finalIteration = F
    while(T) {
      # CLEAN FOLDER IN THE BEGINNING OF EACH ITERATION
      if (!finalIteration)
        do.call(file.remove, list(list.files(paste0(folder_name, sample_name), full.names = TRUE)))
      
      
      
      local_purities <- purities
      local_copy_numbers_used <- copy_numbers_used
      local_cn_states <- cn_states
      local_multipliersDueToLog <- multipliersDueToLog
      local_cnv_states <- statesUsed
      if (sample_name %in% colnames(matrixOfLogFoldOff))
        local_multipliersDueToLogOff <- multipliersDueToLogOff
      
      if (finalIteration ) {
        # if ((abs(max(matrixOfClonality) - min(matrixOfClonality)) > 1000)) {
        #   clonalSignificanceThreshold = 1000
        #   indices = which(matrixOfClonality == min(matrixOfClonality), arr.ind = TRUE)
        #   oneClone = min(matrixOfClonality[indices[1], indices[1]], matrixOfClonality[indices[2], indices[2]])
        #   if (abs(oneClone - matrixOfClonality[indices[1,1], indices[1,2]]) > clonalSignificanceThreshold) {
        #     clonalBestPurities <- rownames(which(matrixOfClonality == min(matrixOfClonality), arr.ind = TRUE))
        #     clonalBestPurities = as.numeric(clonalBestPurities)
        #   } else {
        #     clonalBestPurities <- rownames(which(matrixOfClonality == oneClone, arr.ind = TRUE))
        #     clonalBestPurities = as.numeric(clonalBestPurities)
        #     tableOfBestPurities <- table(rownames(which(matrixOfClonality == oneClone, arr.ind = TRUE)))
        #     clonalBestPurities = as.numeric(names(tableOfBestPurities[which.max(tableOfBestPurities)]))
        #   }
        #   
        #   clonalBestPurities <- c(as.numeric(clonalBestPurities), 0)
        #   indices_to_remove_by_purity <- which(!(purities %in% clonalBestPurities))
        #   local_purities <- purities[-indices_to_remove_by_purity]
        #   local_copy_numbers_used <- copy_numbers_used[-indices_to_remove_by_purity]
        #   local_cn_states <- cn_states[-indices_to_remove_by_purity]
        #   local_multipliersDueToLog <- multipliersDueToLog[-indices_to_remove_by_purity]
        #   local_cnv_states = local_cnv_states[-indices_to_remove_by_purity]
        # }
        clonalBestPurities <- c(as.numeric(clonalBestPurities), 0)
        indices_to_remove_by_purity <- which(!(purities %in% clonalBestPurities))
        local_purities <- purities[-indices_to_remove_by_purity]
        local_copy_numbers_used <- copy_numbers_used[-indices_to_remove_by_purity]
        local_cn_states <- cn_states[-indices_to_remove_by_purity]
        local_multipliersDueToLog <- multipliersDueToLog[-indices_to_remove_by_purity]
        local_cnv_states = local_cnv_states[-indices_to_remove_by_purity]
        if (sample_name %in% colnames(matrixOfLogFoldOff))
          local_multipliersDueToLogOff <- local_multipliersDueToLogOff[-indices_to_remove_by_purity]
      }
      
      # PART FOR MATRIX OF CLONALITY (ONLY 2 CLONES)
      uniqueLocalPurities = unique(local_purities)
      zeroPurity <- which(uniqueLocalPurities == 0)
      uniqueLocalPurities = sort(uniqueLocalPurities[-zeroPurity])
      likeliksFoundCNVsVsPuritiesGlobal = matrix(nrow=0, ncol=length(uniqueLocalPurities))

      
      pvalsForQC <- c()
      threshold = opt$scoreS
      minimum_length_of_CNV = opt$lengthS
      if (!finalIteration) {
        threshold = opt$scoreS + 100
      } else {
        threshold = opt$scoreS
      }
      price_per_tile = 0.1
      initial_state <- 1
      sampleInOfftarget=F
      
      
      localSds = sdsOfProbes * esimtatedVarianceFromSampleNoise[sam_no] * multiplicator
      if (frameworkOff == "offtarget") {
        if (sample_name %in% colnames(matrixOfLogFoldOff)) {
          sam_no_off = which(colnames(matrixOfLogFoldOff) == sample_name)
          localSdsOff = sdsOfProbesOff * esimtatedVarianceFromSampleNoiseOff[sam_no_off] * multiplicatorOff
          sampleInOfftarget = T
        }
      }
      
      
      
      dict_to_output = c()
      
      
      
      
      
      #### CORRECTION - IF THE SAMPLE HAS TOO MANY CNAS, WE EXPECT SOME SHIFT THERE
      if (frameworkDataTypes == "covdepthBAF") {
        sampleName2 <- strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
        position <- which(substring(names(allowedChromsBaf), 1, nchar(sampleName2)) == sampleName2)
        if (length(position) == 1) {
          allowedChromsBafSample <- allowedChromsBaf[[position]]
          
          if (!sampleInOfftarget) {
            allowedChromosomesAutosomesOnly = c()
            for (allowedArm in allowedChromsBafSample) {
              splittedValue <- strsplit(allowedArm, "-")
              chrom = splittedValue[[1]][1]
              if (!chrom %in% c("chrY", "Y", "chrX", "X")) {
                startOfArm = as.numeric(splittedValue[[1]][2])
                endOfArm = as.numeric(splittedValue[[1]][3])
                allowedChromosomesAutosomesOnly = union(allowedChromosomesAutosomesOnly, which(bedFile[,1] == chrom &
                                                                                                 bedFile[,2] >= startOfArm &
                                                                                                 bedFile[,3] <= endOfArm))
              }
            }
            lengthOfRolling = 20
            smoothedLogFold <- rep(0, length(allowedChromosomesAutosomesOnly) - lengthOfRolling)
            matrixOfLogFoldAllowedChrom = matrixOfLogFold[allowedChromosomesAutosomesOnly, sam_no]
            for (i in lengthOfRolling:length(allowedChromosomesAutosomesOnly)) {
              smoothedLogFold[i - lengthOfRolling + 1] = median(matrixOfLogFoldAllowedChrom[(i - lengthOfRolling + 1):i])
            }
            clusteredResult <- densityMclust(smoothedLogFold[which(smoothedLogFold > log2(2/8))])
            print("Mclust finished")
            bigClusters <- which(clusteredResult$parameters$pro > 0.3)
            if (length(bigClusters) == 0) {
              shiftOfCoverage <- median(matrixOfLogFold[allowedChromosomesAutosomesOnly])
            } else {
              shiftOfCoverage = min(clusteredResult$parameters$mean[bigClusters])
            }
          } else {
            sam_no_off = which(colnames(matrixOfLogFoldOff) == sample_name)
            globalBed <- rbind(bedFile, bedFileOfftarget)
            globalLogFold <- c( matrixOfLogFold[,sam_no], matrixOfLogFoldOff[,sam_no_off])
            allowedChromosomesAutosomesOnly = c()
            for (allowedArm in allowedChromsBafSample) {
              splittedValue <- strsplit(allowedArm, "-")
              chrom = splittedValue[[1]][1]
              if (!chrom %in% c("chrY", "Y", "chrX", "X")) {
                startOfArm = as.numeric(splittedValue[[1]][2])
                endOfArm = as.numeric(splittedValue[[1]][3])
                allowedChromosomesAutosomesOnly = union(allowedChromosomesAutosomesOnly, which(globalBed[,1] == chrom &
                                                                                                 globalBed[,2] >= startOfArm &
                                                                                                 globalBed[,3] <= endOfArm))
              }
            }
            lengthOfRolling = 20
            smoothedLogFold <- rep(0, length(allowedChromosomesAutosomesOnly) - lengthOfRolling)
            globalLogFoldAllowedChroms = globalLogFold[allowedChromosomesAutosomesOnly]
            for (i in lengthOfRolling:length(allowedChromosomesAutosomesOnly)) {
              smoothedLogFold[i - lengthOfRolling + 1] = median(globalLogFoldAllowedChroms[(i - lengthOfRolling + 1):i])
            }
            clusteredResult <- densityMclust(smoothedLogFold[which(smoothedLogFold > log2(2/8))])
            print("Mclust finished")
            bigClusters <- which(clusteredResult$parameters$pro > 0.3)
            if (length(bigClusters) == 0) {
              shiftOfCoverage <- median(globalLogFold[allowedChromosomesAutosomesOnly])
            } else {
              shiftOfCoverage = min(clusteredResult$parameters$mean[bigClusters])
            }
          }
          matrixOfLogFold[,sam_no] = matrixOfLogFold[,sam_no] - shiftOfCoverage
          if (sampleInOfftarget)
            matrixOfLogFoldOff[,sam_no_off] = matrixOfLogFoldOff[,sam_no_off] - shiftOfCoverage
        }
        
      }
      
      
      matrixOfLogFoldCorrectedSmall = matrixOfLogFold
      matrixOfLogFoldCorrectedSmall[which(matrixOfLogFoldCorrectedSmall < log2(min(local_cn_states) / 2))] = log2(min(local_cn_states) / 2)
      matrixOfLogFoldCorrectedSmall[which(matrixOfLogFoldCorrectedSmall > log2(max(local_cn_states) / 2))] = log2(max(local_cn_states) / 2)
      
      matrix_of_likeliks <- form_matrix_of_likeliks_one_sample(1, ncol(matrixOfLogFoldCorrectedSmall), matrixOfLogFoldCorrectedSmall[,sam_no], localSds, log2(local_cn_states/2), local_multipliersDueToLog)
      
      if (genderOfSamples[germline_sample_no] == "M") {
        matrix_of_likeliks[which(bedFile[,1] %in% c("chrX","chrY")),] = form_matrix_of_likeliks_one_sample(
          1, ncol(matrixOfLogFoldCorrectedSmall), matrixOfLogFoldCorrectedSmall[which(bedFile[,1] %in% c("chrX","chrY")),sam_no], 
          localSds[which(bedFile[,1] %in% c("chrX","chrY"))], log2((1 - local_purities) + local_purities * local_copy_numbers_used), local_multipliersDueToLog)
      }
      
      matrOfSNVlikeliks <- matrix(0, nrow=0, ncol=length(local_purities))
      
      
      
      ### ADD LIKELIHOODS
      if (frameworkDataTypes == "covdepthBAF") {
        print("Started BAF calculation")
        print(Sys.time())
        
        numberOfAssignedPositions = 0
        sampleName2 <- strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
        position <- which(substring(names(allowedChromsBaf), 1, nchar(sampleName2)) == sampleName2)
        if (length(position) == 1) {
          bAlleleFreqsTumor <- bAlleleFreqsAllSamples[[position]][[ strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1] ]]
          bAlleleFreqsNormal <- bAlleleFreqsAllSamples[[position]][[ strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2] ]]
          
          # calculate median correction factor
          allowedChromosomesAutosomesOnly = which(!bAlleleFreqsTumor[,1] %in% c("X","Y","chrX","chrY"))
          multiplierOfSNVsDueToMapping <- median(as.numeric(bAlleleFreqsNormal[allowedChromosomesAutosomesOnly,5]))
          print("Multiplier of allele balance of a particular sample")
          print(multiplierOfSNVsDueToMapping)
          
          numOfSNVs = nrow(bAlleleFreqsTumor)
          reduceOfSNVsSize = 1
          if (numOfSNVs > 2000) {
            for (reduceOfSNVsSize in 2:100) {
              if (numOfSNVs / reduceOfSNVsSize < 2000 ) {
                break
              }
            }
            reduceOfSNVsSize = reduceOfSNVsSize - 1
          }
          if (length(closestBedRegions) == 0) closestBedRegions = rep(0, nrow(bAlleleFreqsTumor))
          for (i in 1:nrow(bAlleleFreqsTumor)) {
            # To avoid computationally expensive steps on the start of estimation
            if (i %% 100 == 0) {
              print(i / nrow(bAlleleFreqsTumor))
              print(Sys.time())
            }
            if (vectorsWithRegionCoordsFilled) {
              closestBedRegion = closestBedRegions[i]
            } else {
              closestBedRegion <- which(bedFile[,1] == bAlleleFreqsTumor[i,1] & bedFile[,2] <= bAlleleFreqsTumor[i,2] & bedFile[,3] >= bAlleleFreqsTumor[i,3])
              if (length(closestBedRegion) >= 1) {
                closestBedRegion = closestBedRegion[1]
                closestBedRegions[i] = closestBedRegion
              } else {
                # if there is no matching position we put 0
                closestBedRegions[i] = 0
                closestBedRegion = 0
              }
            }
            if (!finalIteration) {
              if (i %% reduceOfSNVsSize != 0 & reduceOfSNVsSize != 1) next
            }
            altAlleleDepth <- as.numeric(bAlleleFreqsTumor[i,5])
            overallDepth <- round(as.numeric(bAlleleFreqsTumor[i,6]))
            altAlleleDepth = round(altAlleleDepth * overallDepth)
            
            if (length(closestBedRegion) == 1 & closestBedRegion != 0) {
              numberOfAssignedPositions = numberOfAssignedPositions + 1
              
              pList = list()
              vecOfLikeliks <- rep(0, ncol(matrix_of_likeliks))
              for (j in 1:length(local_cn_states)) {
                pur = local_purities[j]
                cn = local_copy_numbers_used[j]
                stateUsed = local_cnv_states[j]
                listOfLikelikAndPList = likelihoodOfSNVBasedOnCN(altAlleleDepth, overallDepth, pur, cn, stateUsed, multiplierOfSNVsDueToMapping, pList)
                
                likelihood = -2 * listOfLikelikAndPList[[1]]
                if (listOfLikelikAndPList[[2]])
                  pList = listOfLikelikAndPList[[3]]
                vecOfLikeliks[j] = likelihood
              }
              
              oldLikeliks <- matrix_of_likeliks[i,] 
              matrix_of_likeliks[closestBedRegion,] = matrix_of_likeliks[closestBedRegion,] + vecOfLikeliks
            }
          }
          vectorsWithRegionCoordsFilled = T
        }
        print("Finished BAF calculation")
        print(Sys.time())
      }
      # Reduce probability of unrealistic state = only super strong evidence is required to go for them
      fine_for_unrealistic_state = 0.5
      set_of_unrealistic_states = c("LOHDup", "CNVboth", "CNVcomplex")
      whichAreUnrealistic <- which(local_cnv_states %in% set_of_unrealistic_states)
      matrix_of_likeliks[,whichAreUnrealistic] = matrix_of_likeliks[,whichAreUnrealistic] + fine_for_unrealistic_state
      
      sizesOfPointsFromLocalSds <- 0.5 / localSds 
      if (sampleInOfftarget) {
        matrixOfLogFoldOffCorrectedExtraSmallValues <- matrixOfLogFoldOff
        matrixOfLogFoldOffCorrectedExtraSmallValues[which(matrixOfLogFoldOffCorrectedExtraSmallValues < log2(min(local_cn_states) / 2))] = log2(min(local_cn_states) / 2)
        matrixOfLogFoldOffCorrectedExtraSmallValues[which(matrixOfLogFoldOffCorrectedExtraSmallValues > log2(max(local_cn_states) / 2))] = log2(max(local_cn_states) / 2)
        matrix_of_likeliks_off <- form_matrix_of_likeliks_one_sample(1, ncol(matrixOfLogFoldOffCorrectedExtraSmallValues), matrixOfLogFoldOffCorrectedExtraSmallValues[,sam_no_off], localSdsOff, log2(local_cn_states/2), local_multipliersDueToLogOff)
        if (genderOfSamples[germline_sample_no] == "M") {
          matrix_of_likeliks_off[which(bedFileOfftarget[,1] %in% c("chrX","chrY")),] = form_matrix_of_likeliks_one_sample(
            1, ncol(matrixOfLogFoldOffCorrectedExtraSmallValues), matrixOfLogFoldOffCorrectedExtraSmallValues[which(bedFileOfftarget[,1] %in% c("chrX","chrY")),sam_no_off], 
            localSdsOff[which(bedFileOfftarget[,1] %in% c("chrX","chrY"))], log2((1 - local_purities) + local_purities * local_copy_numbers_used), local_multipliersDueToLog)
        }
        
        
        
        globalMatrOfLikeliks <- rbind(matrix_of_likeliks, matrix_of_likeliks_off)
        globalBed <- rbind(bedFile, bedFileOfftarget)
        sizesOfPointsFromLocalSdsOff <- 0.5 / localSdsOff
        vecOfOrder = order(globalBed[,1], as.numeric(globalBed[,2]))
        globalSizesOfPoints <- c(sizesOfPointsFromLocalSds, sizesOfPointsFromLocalSdsOff)[vecOfOrder]
        globalMatrOfLikeliks <- globalMatrOfLikeliks[vecOfOrder,]
        globalLogFold <- c( matrixOfLogFold[,sam_no], matrixOfLogFoldOff[,sam_no_off])[vecOfOrder]
        globalSds <-  c(localSds, localSdsOff)[vecOfOrder]
        globalBed <- globalBed[vecOfOrder,]
      }
      
      
      
      
      
      found_CNVs_total <- matrix(0, nrow=0, ncol=10)
      colnames(found_CNVs_total) <- c("#chr", "start", "end", "tumor_CN_change", "tumor_clonality", "CN_change", "loglikelihood", "number_of_regions", "state", "genes")
      allDetectedPurities = c()
      for (l in 1:length(left_borders)) {
        
        
        
        
        chrom = names(left_borders)[l]
        
        
        if (chrom == "chrY" & genderOfSamples[sam_no_germline] == "F") {
          next
        } else if (genderOfSamples[sam_no_germline] == "M") {
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
              which_to_allow = which(globalBed[,1] == chrom & globalBed[,2] <= start )
            } else {
              which_to_allow = which(globalBed[,1] == chrom & globalBed[,2] >= end )
            }
            toyBedFile = globalBed[which_to_allow,]
            
            toyMatrixOfLikeliks = globalMatrOfLikeliks[which_to_allow,]
            toyLogFoldChange = globalLogFold[which_to_allow]
            
            found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, 1))
            
            toySds <- globalSds[which_to_allow]
            
            #if (pvalueForThisArmQC > 0) {
            #pvalsForQC <- c(pvalsForQC, pvalueForThisArmQC)
            #}
            
            
            toySizesOfPointsFromLocalSds = c(sizesOfPointsFromLocalSds, sizesOfPointsFromLocalSdsOff)[vecOfOrder]
            toySizesOfPointsFromLocalSds = toySizesOfPointsFromLocalSds[which_to_allow]
          } else {
            if (k == 1) {
              which_to_allow = which(bedFile[,1] == chrom & bedFile[,2] <= start )
            } else {
              which_to_allow = which(bedFile[,1] == chrom & bedFile[,2] >= end )
            }
            toyMatrixOfLikeliks = matrix_of_likeliks[which_to_allow,]
            toyBedFile = bedFile[which_to_allow,]
            found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, 1))
            toySds <- localSds[which_to_allow]
            
            
            # if (pvalueForThisArmQC > 0) {
            #pvalsForQC <- c(pvalsForQC, pvalueForThisArmQC)
            #}
            
            toyLogFoldChange = matrixOfLogFold[which_to_allow, sam_no]
            toySizesOfPointsFromLocalSds = sizesOfPointsFromLocalSds[which_to_allow]
          }
          
          if (nrow(found_CNVs) > 0 & !chrom %in% c("chrX", "chrY", "X", "Y")) {
            likeliksFoundCNVsVsPurities <- matrix(nrow=nrow(found_CNVs), ncol=length(uniqueLocalPurities))
            for (m in 1:length(uniqueLocalPurities)) {
              localPurityCurrent = uniqueLocalPurities[m]
              for (q in 1:nrow(found_CNVs)) {
                startOfCNV <- found_CNVs[q,2]
                endOfCNV <- found_CNVs[q,3]
                if (endOfCNV - startOfCNV > 3) { 
                  penalty = penaltyForHigherCN * abs(2 - local_copy_numbers_used[which(local_purities == localPurityCurrent)])
                  likeliksFoundCNVsVsPurities[q, m] = min(apply(toyMatrixOfLikeliks[(startOfCNV + 1):(endOfCNV - 1),which(local_purities == localPurityCurrent)], 2, sum) + penalty)
                }
              }
            }
            likeliksFoundCNVsVsPuritiesGlobal = rbind(likeliksFoundCNVsVsPuritiesGlobal, likeliksFoundCNVsVsPurities)
          }
          
          
          
          ### IGV PLOTTING
          if(opt$debug) {
            print("START OF IGV PLOTTING")
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
                                                outputFileNameDots, sample_name, toyLogFoldChange, reverseFunctionUsedToTransform, local_cn_states)
          }
          if(opt$debug) {
            print("END OF IGV PLOTTING")
          }
          ### END OF IGV PLOTTING
          
          
          
          if (nrow(found_CNVs) == 0 & length(which_to_allow) > 1) {
            found_CNVs = matrix(c(-1000, 1, length(which_to_allow), 1), nrow=1)
          }
          
          
          
          if (nrow(found_CNVs) > 0) {
            cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyLogFoldChange, toyBedFile, output_of_plots, chrom, local_cn_states, local_copy_numbers_used, local_purities, local_cnv_states, 
                                            toySizesOfPointsFromLocalSds, plottingOfPNGs)
            if (found_CNVs[1,1] != -1000) {
              found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
            }
          }
        }
        
      }
      
      
      if (finalIteration == T) break
      
      
      
      
      
      
      if (nrow(likeliksFoundCNVsVsPuritiesGlobal) > 0) {
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
        minPointToNormalise = quantile(matrixOfClonality, 0.75)
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
        maxNumOfClones <- 7
        localPurityStates = 1:length(uniqueLocalPurities)
        resultBestCombination = 0
        for (m in 1:maxNumOfClones) {
          combinationsOfPurities <- combn(localPurityStates, m)
          bestCombination = 1
          minResult = 0
          for (r in 1:nrow(likeliksFoundCNVsVsPuritiesGlobal)) {
            minResult = minResult + min(likeliksFoundCNVsVsPuritiesGlobal[r,combinationsOfPurities[,1]])
          }
          
          for (q in 1:ncol(combinationsOfPurities)) {
            minResultForCombination = 0
            for (r in 1:nrow(likeliksFoundCNVsVsPuritiesGlobal)) {
              minResultForCombination = minResultForCombination + min(
              likeliksFoundCNVsVsPuritiesGlobal[r,combinationsOfPurities[,q]])
            }
            if (minResult > minResultForCombination) {
              bestCombination = q
              minResult = minResultForCombination
            }
          }
          if (m == 1) {
            minResultSoFar = minResult
            resultBestCombination = combinationsOfPurities[bestCombination]
          } else {
            if (minResult > minResultSoFar - opt$scoreS - 100) {
              print(uniqueLocalPurities[resultBestCombination])
              break
            } else {
              print("Current improve")
              print(minResultSoFar - minResult)
              minResultSoFar = minResult
              resultBestCombination = combinationsOfPurities[,bestCombination]
            }
          }
        }
        clonalBestPurities = uniqueLocalPurities[resultBestCombination]
        if (length(clonalBestPurities) == 0) {
          clonalBestPurities = c(0, 1)
        }
        
      } else {
        print("No CNVs found in this sample for finding clonality.")
        clonalBestPurities = c(0, 1)
      }
      finalIteration = T
    }
    
    if (length(pvalsForQC > 1)) {
      finalPValue <- 0
    } else {
      finalPValue = 0
    }
    fileToOut <- paste0(folder_name, sample_name, paste0("/CNAs_", sample_name, ".txt"))
    fileConn<-file(fileToOut)
    writeLines(c(paste("##"," QC ", 0, "clonality by BAF (if != 1):", paste(round(unique(local_purities), digits=3), collapse=";"), collapse = " ")), fileConn)
    close(fileConn)
    if(opt$debug) {
      print(found_CNVs_total)
    }
    write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)	
  }
}


