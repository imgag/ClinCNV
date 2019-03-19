
no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
cl<-makeCluster(no_cores, type="FORK")
registerDoParallel(cl)


### PROCESSING OF SOMATIC VARIANTS
setwd(opt$folderWithScript)
source(paste0(opt$folderWithScript, "/somatic/helpersSomatic.R"),local=TRUE)

vect_of_norm_likeliks <- fast_dt_list(as.numeric(opt$degreesOfFreedomStudent))
vect_of_t_likeliks <- fast_dt_list(as.numeric(opt$degreesOfFreedomStudent))


print(paste("Work on data preparation for somatic samples started (log-fold change matrices plus parameters estimation)", Sys.time()))

listOfValue <- formilngLogFoldChange(pairs, tmpNormal, tumor, bedFileForCluster, genderOfSamples)
matrixOfLogFold <- listOfValue[[1]]
gendersInOntargetMatrix = listOfValue[[2]]



matrixWithSdsList <- findSDsOfSamples(pairs, tmpNormal, tumor, bedFileForCluster, bordersOfChroms, genderOfSamples)
matrixWithSds = matrixWithSdsList[[1]]
sdsOfSomaticSamples <- matrixWithSds[4,]
sdsOfProbes <- matrixWithSdsList[[2]]

# QC control
sampleLevelOfNoise = apply(matrixOfLogFold[which(!bedFileForCluster[,1] %in% c("chrX","chrY")),], 2, Sn)
samplesNotPassedQC = colnames(matrixOfLogFold)[which(sampleLevelOfNoise > 1.0)]
if (length(samplesNotPassedQC) > 0) {
  print("Samples don't pass our QC - too noisy!")
  print(samplesNotPassedQC)
  sdsOfSomaticSamples = sdsOfSomaticSamples[-which(colnames(matrixOfLogFold) %in% samplesNotPassedQC)]
  matrixWithSds = matrixWithSds[,-which(colnames(matrixOfLogFold) %in% samplesNotPassedQC)]
  gendersInOntargetMatrix = gendersInOntargetMatrix[-which(colnames(matrixOfLogFold) %in% samplesNotPassedQC)]
  matrixOfLogFold = matrixOfLogFold[,-which(colnames(matrixOfLogFold) %in% samplesNotPassedQC)]
}

## QC FOR PROBES
probesToRemove <- probeLevelQC(matrixOfLogFold, sdsOfProbes, sdsOfSomaticSamples, gendersInOntargetMatrix, bedFileForCluster)
if (length(probesToRemove) > 0) {
  sdsOfProbes = sdsOfProbes[-probesToRemove]
  matrixOfLogFold = matrixOfLogFold[-probesToRemove,]
  bedFileForCluster = bedFileForCluster[-probesToRemove,]
}

if (frameworkOff == "offtarget") {
  listOfValueOff <- formilngLogFoldChange(pairs, normalOff[,which(colnames(normalOff) %in% colnames(tmpNormal))], tumorOff, bedFileForClusterOff, genderOfSamples)
  matrixOfLogFoldOff =  listOfValueOff[[1]]
  gendersInOfftargetMatrix = listOfValueOff[[2]]

  if (length(samplesNotPassedQC) > 0) {
    if (length(which(colnames(matrixOfLogFoldOff) %in% samplesNotPassedQC)) > 0) {
      gendersInOfftargetMatrix = gendersInOfftargetMatrix[-which(colnames(matrixOfLogFoldOff) %in% samplesNotPassedQC)]
      matrixOfLogFoldOff = matrixOfLogFoldOff[,-which(colnames(matrixOfLogFoldOff) %in% samplesNotPassedQC)]
    }
  }
  

  
  
  bordersOfChroms <- getBordersOfChromosomes(bedFileForClusterOff)
  matrixWithSdsOffList <- findSDsOfSamples(pairs, normalOff[,which(colnames(normalOff) %in% colnames(tmpNormal))], tumorOff, bedFileForClusterOff, bordersOfChroms, genderOfSamples)
  matrixWithSdsOff = matrixWithSdsOffList[[1]]
  sdsOfSomaticSamplesOff <- matrixWithSdsOff[4,]
  sdsOfProbesOff <- matrixWithSdsOffList[[2]]
  
  probesToRemove <- probeLevelQC(matrixOfLogFoldOff, sdsOfProbesOff, sdsOfSomaticSamplesOff, gendersInOfftargetMatrix, bedFileForClusterOff)
  if (length(probesToRemove) > 0) {
    sdsOfProbesOff = sdsOfProbesOff[-probesToRemove]
    matrixOfLogFoldOff = matrixOfLogFoldOff[-probesToRemove,]
    bedFileForClusterOff = bedFileForClusterOff[-probesToRemove,]
  }
  #sdsOfProbesOff <- sapply(1:nrow(matrixOfLogFoldOff), function(i) {determineSDsOfSomaticProbe(matrixOfLogFoldOff[i,], i)})
  
  #listOfVarianceAndMultiplicatorOff <- esimtateVarianceFromSampleNoise(sdsOfSomaticSamplesOff, 10000)
  #esimtatedVarianceFromSampleNoiseOff <- listOfVarianceAndMultiplicatorOff[[2]]
  #multiplicatorOff <- listOfVarianceAndMultiplicatorOff[[1]]
  
}






cn_states <- c()
copy_numbers = 0:15
purity <- seq(from=opt$minimumPurity, to=100.1, by=opt$purityStep) / 100
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
      statesUsed <- c(statesUsed, "CNVcomplex2")
    }
    if (cn >= 7 & cn <= 8) {
      cn_state <- (1 - pur) * 2 + pur * cn
      cn_states <- c(cn_states, cn_state)
      copy_numbers_used <- c(copy_numbers_used, cn)
      purities <- c(purities, pur)
      statesUsed <- c(statesUsed, "CNVcomplex3")
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







folder_name <- paste0(opt$out, "/somatic/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

allPotentialPurities <- unique(purities)
penaltyForHigherCN = 10
clonalityForChecking = 0.4
print(paste("Work on actual calling started.", Sys.time()))

normalNames = sapply(1:length(allowedChromsBaf), function(i) {strsplit(names(allowedChromsBaf)[i], split="-")[[1]][2]})
for (sam_no in 1:ncol(matrixOfLogFold)) {
  sample_name <- colnames(matrixOfLogFold)[sam_no]

  germline_sample_no = which(colnames(tmpNormal) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2])
  tumor_sample_no = which(colnames(tumor) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1])
  if (strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1] %in% colnames(tumorOff)) {
    tumor_sample_no_off = which(colnames(tumorOff) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1])
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
  
  if (!dir.exists(paste0(folder_name, sample_name)) | (!is.null(opt$reanalyseCohort))) {
    
    dir.create(paste0(folder_name, sample_name))
    setwd(paste0(folder_name, sample_name))
    
    copyNumbersInsideExpectedPurities = F
    
    emptyChroms = matrix(0, ncol=2, nrow=0)
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
        threshold = opt$scoreS
      } else {
        threshold = opt$scoreS
      }
      price_per_tile = 0.1
      initial_state <- 1
      sampleInOfftarget=F
      
      
      localSds = sdsOfProbes * (sdsOfSomaticSamples[sam_no])
      if (frameworkOff == "offtarget") {
        if (sample_name %in% colnames(matrixOfLogFoldOff)) {
          sam_no_off = which(colnames(matrixOfLogFoldOff) == sample_name)
          localSdsOff = sdsOfProbesOff * (sdsOfSomaticSamplesOff[sam_no_off])
          sampleInOfftarget = T
        }
      }
      
      
      
      dict_to_output = c()
      
      
      
      
      
      #### CORRECTION - IF THE SAMPLE HAS TOO MANY CNAS, WE EXPECT SOME SHIFT THERE
      if (frameworkDataTypes == "covdepthBAF" & sample_name %in% normalNames) {
        sampleName2 <- strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
        tumorNames = sapply(1:length(allowedChromsBaf), function(i) {strsplit(names(allowedChromsBaf)[i], split="-")[[1]][1]})
        position <- which(tumorNames == sampleName2)
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
                allowedChromosomesAutosomesOnly = union(allowedChromosomesAutosomesOnly, which(bedFileForCluster[,1] == chrom &
                                                                                                 bedFileForCluster[,2] >= startOfArm &
                                                                                                 bedFileForCluster[,3] <= endOfArm))
              }
            }
            lengthOfRolling = 30
            matrixOfLogFoldAllowedChrom = matrixOfLogFold[allowedChromosomesAutosomesOnly, sam_no]
            
            smoothedLogFold = runmed(matrixOfLogFoldAllowedChrom, k = lengthOfRolling)
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
            globalBed <- rbind(bedFileForCluster, bedFileForClusterOff)
            globalLogFold <- c( matrixOfLogFold[,sam_no], matrixOfLogFoldOff[,sam_no_off])
            allowedChromosomesAutosomesOnly = c()
            smoothedLogFold= c()
            for (allowedArm in allowedChromsBafSample) {
              splittedValue <- strsplit(allowedArm, "-")
              chrom = splittedValue[[1]][1]
              if (!chrom %in% c("chrY", "Y", "chrX", "X")) {
                startOfArm = as.numeric(splittedValue[[1]][2])
                endOfArm = as.numeric(splittedValue[[1]][3])
                lengthOfRolling = min(100, round((endOfArm - startOfArm)/5))
                allowedChromosomesAutosomesOnly = union(allowedChromosomesAutosomesOnly, which(globalBed[,1] == chrom &
                                                                                                 globalBed[,2] >= startOfArm &
                                                                                                 globalBed[,3] <= endOfArm))
                smoothedLogFold = c(smoothedLogFold, runmed(globalLogFold[which(globalBed[,1] == chrom &
                                                                                   globalBed[,2] >= startOfArm &
                                                                                   globalBed[,3] <= endOfArm)], k = lengthOfRolling))
              }
            }
            #globalLogFoldAllowedChroms = globalLogFold[allowedChromosomesAutosomesOnly]
            #smoothedLogFold = runmed(globalLogFoldAllowedChroms, k = lengthOfRolling)
            clusteredResult <- densityMclust(smoothedLogFold[which(smoothedLogFold > log2(3/8))], model="E")
            print("Mclust finished")
            bigClusters <- which(clusteredResult$parameters$pro > 0.25)
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
        if (length(which(bedFileForCluster[,1] %in% c("chrX","chrY"))) > 0)
          matrix_of_likeliks[which(bedFileForCluster[,1] %in% c("chrX","chrY")),] = form_matrix_of_likeliks_one_sample(
            1, ncol(matrixOfLogFoldCorrectedSmall), matrixOfLogFoldCorrectedSmall[which(bedFileForCluster[,1] %in% c("chrX","chrY")),sam_no], 
            localSds[which(bedFileForCluster[,1] %in% c("chrX","chrY"))], log2((1 - local_purities) + local_purities * local_copy_numbers_used), local_multipliersDueToLog)
      }
      
      
      ### ADD LIKELIHOODS
      bAlleleFreqsTumor = NULL
      bAlleleFreqsNormal = NULL
      if (frameworkDataTypes == "covdepthBAF" & sample_name %in% normalNames) {
        print("Started BAF calculation")
        print(Sys.time())
        
        numberOfAssignedPositions = 0
        sampleName2 <- strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
        tumorNames = sapply(1:length(allowedChromsBaf), function(i) {strsplit(names(allowedChromsBaf)[i], split="-")[[1]][1]})
        position <- which(tumorNames == sampleName2)
        if (length(position) == 1) {
          bAlleleFreqsTumor <- bAlleleFreqsAllSamples[[position]][[ strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1] ]]
          bAlleleFreqsNormal <- bAlleleFreqsAllSamples[[position]][[ strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2] ]]
          if (genderOfSamples[germline_sample_no] == "M") {
            bAlleleFreqsTumor = bAlleleFreqsTumor[which(!bAlleleFreqsTumor[,1] %in% c("chrX", "chrY")),]
            bAlleleFreqsNormal = bAlleleFreqsNormal[which(!bAlleleFreqsNormal[,1] %in% c("chrX", "chrY")),]
          }
          overdispersionNormal = overdispersionsNormal[[position]]
          overdispersionTumor = overdispersionsTumor[[position]]
          # calculate median correction factor
          allowedChromosomesAutosomesOnly = which(!bAlleleFreqsTumor[,1] %in% c("X","Y","chrX","chrY"))
          multiplierOfSNVsDueToMapping <- median(as.numeric(bAlleleFreqsNormal[allowedChromosomesAutosomesOnly,5]))
          print("Multiplier of allele balance of a particular sample")
          print(multiplierOfSNVsDueToMapping)
          
          numOfSNVs = nrow(bAlleleFreqsTumor)
          if (length(closestBedRegions) == 0) closestBedRegions = rep(0, nrow(bAlleleFreqsTumor))
          for (i in 1:nrow(bAlleleFreqsTumor)) {
            # To avoid computationally expensive steps on the start of estimation
            if (vectorsWithRegionCoordsFilled) {
              closestBedRegion = closestBedRegions[i]
            } else {
              closestBedRegion <- which(bedFileForCluster[,1] == bAlleleFreqsTumor[i,1] & bedFileForCluster[,2] <= bAlleleFreqsTumor[i,2] & bedFileForCluster[,3] >= bAlleleFreqsTumor[i,3])
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
            coordsIncludedAtFirst = setdiff(1:nrow(bAlleleFreqsTumor), union(allowedChromosomesAutosomesOnlySNV, which(bAlleleFreqsTumor[,1] %in% ifelse(genderOfSamples[germline_sample_no] == "M",
                                                                                                                                                         c("chrY"),
                                                                                                                                                         c("chrX", "chrY")))))
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
          matrixOfBAFLikeliks = foreach (i = 1:nrow(bAlleleFreqsTumorToy), .combine='rbind') %dopar% {
            altAlleleDepth <- as.numeric(bAlleleFreqsTumorToy[i,5])
            overallDepth <- round(as.numeric(bAlleleFreqsTumorToy[i,6]))
            altAlleleDepth = round(altAlleleDepth * overallDepth)
            overdispersionValue = overdispersionTumorToy[i]
            
            if (closestBedRegionsToy[i] != 0) {
              numberOfAssignedPositions = numberOfAssignedPositions + 1
              
              pList = list()
              vecOfLikeliks <- rep(0, ncol(matrix_of_likeliks))
              for (j in 1:length(local_cn_states)) {
                pur = local_purities[j]
                cn = local_copy_numbers_used[j]
                stateUsed = local_cnv_states[j]
                
                listOfLikelikAndPList = likelihoodOfSNVBasedOnCN(altAlleleDepth, overallDepth, pur, cn, stateUsed, multiplierOfSNVsDueToMapping, pList, overdispersionValue)
                
                likelihood = -2 * listOfLikelikAndPList[[1]]
                if (listOfLikelikAndPList[[2]])
                  pList = listOfLikelikAndPList[[3]]
                vecOfLikeliks[j] = likelihood
              }
              vecOfLikeliks
            } else {
              return(rep(0, length(local_cn_states)))
            }
          }
          for (i in 1:nrow(bAlleleFreqsTumor)) {
            closestBedRegion = closestBedRegionsToy[i]
            if (!is.na(closestBedRegion))
              if (closestBedRegion != 0)
                matrix_of_likeliks[closestBedRegion,] = matrix_of_likeliks[closestBedRegion,] + matrixOfBAFLikeliks[i,]
          }
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
        matrix_of_likeliks_off <- form_matrix_of_likeliks_one_sample(1, ncol(matrixOfLogFoldOffCorrectedExtraSmallValues), matrixOfLogFoldOffCorrectedExtraSmallValues[,sam_no_off], localSdsOff, log2(local_cn_states/2), local_multipliersDueToLogOff)
        if (genderOfSamples[germline_sample_no] == "M") {
          matrix_of_likeliks_off[which(bedFileForClusterOff[,1] %in% c("chrX","chrY")),] = form_matrix_of_likeliks_one_sample(
            1, ncol(matrixOfLogFoldOffCorrectedExtraSmallValues), matrixOfLogFoldOffCorrectedExtraSmallValues[which(bedFileForClusterOff[,1] %in% c("chrX","chrY")),sam_no_off], 
            localSdsOff[which(bedFileForClusterOff[,1] %in% c("chrX","chrY"))], log2((1 - local_purities) + local_purities * local_copy_numbers_used), local_multipliersDueToLog)
        }
        
        
        
        globalMatrOfLikeliks <- rbind(matrix_of_likeliks, matrix_of_likeliks_off)
        globalBed <- rbind(bedFileForCluster, bedFileForClusterOff)
        sizesOfPointsFromLocalSdsOff <- 0.5 / localSdsOff
        vecOfOrder = order(globalBed[,1], as.numeric(globalBed[,2]))
        globalSizesOfPoints <- c(sizesOfPointsFromLocalSds, sizesOfPointsFromLocalSdsOff)[vecOfOrder]
        globalMatrOfLikeliks <- globalMatrOfLikeliks[vecOfOrder,]
        globalLogFold <- c( matrixOfLogFold[,sam_no], matrixOfLogFoldOff[,sam_no_off])[vecOfOrder]
        globalSds <-  c(localSds, localSdsOff)[vecOfOrder]
        globalBed <- globalBed[vecOfOrder,]
      }
      print(paste("Block before CNV detection finished", Sys.time()))
      
      
      
      
      found_CNVs_total <- matrix(0, nrow=0, ncol=10)
      colnames(found_CNVs_total) <- c("#chr", "start", "end", "tumor_CN_change", "tumor_clonality", "CN_change", "loglikelihood", "number_of_regions", "state", "genes")
      allDetectedPurities = c()
      for (l in 1:length(left_borders)) {
        
        
        
        
        chrom = names(left_borders)[l]
        print(paste(chrom, Sys.time()))
        
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
            } else {
              which_to_allow = which(globalBed[,1] == chrom & as.numeric(globalBed[,2]) >= as.numeric(end) )
            }
            toyBedFile = globalBed[which_to_allow,]
            
            toyMatrixOfLikeliks = globalMatrOfLikeliks[which_to_allow,]
            toyLogFoldChange = globalLogFold[which_to_allow]
            
            toySds <- globalSds[which_to_allow]
            
            toySizesOfPointsFromLocalSds = c(sizesOfPointsFromLocalSds, sizesOfPointsFromLocalSdsOff)[vecOfOrder]
            toySizesOfPointsFromLocalSds = toySizesOfPointsFromLocalSds[which_to_allow]
          } else {
            if (k == 1) {
              which_to_allow = which(bedFileForCluster[,1] == chrom & bedFileForCluster[,2] <= start )
            } else {
              which_to_allow = which(bedFileForCluster[,1] == chrom & bedFileForCluster[,2] >= end )
            }
            toyMatrixOfLikeliks = matrix_of_likeliks[which_to_allow,]
            toyBedFile = bedFileForCluster[which_to_allow,]
            toyLogFoldChange = matrixOfLogFold[which_to_allow, sam_no]
            

            toySds <- localSds[which_to_allow]
            
            toySizesOfPointsFromLocalSds = sizesOfPointsFromLocalSds[which_to_allow]
          }
          blocked_states = c()
          if (length(toyLogFoldChange) > opt$lengthS) {
            arrayOfMediansOfToyLogFold = runmed(toyLogFoldChange, round(opt$lengthS/ 2))
            #arrayOfMediansOfToyLogFold <- sapply(1:(length(toyLogFoldChange) - opt$lengthS), function(i) {median(toyLogFoldChange[i:(i + opt$lengthS)])})
            if (!finalIteration) {
              diffsFromCoverage <- sapply(1:length(local_cn_states), function(i) {min(abs(log2(local_cn_states[i] / local_cn_states[initial_state]) - (arrayOfMediansOfToyLogFold)))})
              blocked_states = c(setdiff(c(1,2), initial_state),
                                 which(diffsFromCoverage > 0.1))
            } else {
            blocked_states = c(setdiff(c(1,2), initial_state),
                               which(log2(local_cn_states / local_cn_states[initial_state]) < min(arrayOfMediansOfToyLogFold) - 0.1 | log2(local_cn_states / local_cn_states[initial_state]) > max(arrayOfMediansOfToyLogFold) + 0.1))
            }
            if (initial_state %in% blocked_states) {
              blocked_states = blocked_states[-which(blocked_states == initial_state)]
            }
            if (length(local_cn_states) - length(blocked_states) == 1) {
              blocked_states <- c()
            }
          }
          # BLOCK WITH PENALTIES
          copy_numbers_for_penalties = 3 - local_copy_numbers_used
          copy_numbers_for_penalties[which(copy_numbers_for_penalties > 0)] = 0
          set_of_unrealistic_states = c("LOHDup", "CNVcomplex2", "CNVcomplex3")
          whichAreUnrealistic <- which(local_cnv_states %in% set_of_unrealistic_states)
          penalties = penaltyForHigherCN * abs(copy_numbers_for_penalties)
          penalties[whichAreUnrealistic] = penalties[whichAreUnrealistic] + 20
          if (length(local_cn_states) - length(blocked_states) > 1) {
            found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, 1, blocked_states, penalties))
          } else {
            found_CNVs = matrix(0, nrow=0, ncol=10)
          }
          
          # BAFs from this chromosome
          if (frameworkDataTypes == "covdepthBAF" & !is.null(overdispersionNormal) & nrow(found_CNVs) > 0 & sampe_name %in% normalNames) {
            bafsFromThisChr = which(bAlleleFreqsNormal[,1] == chrom)
            listOfCNVsThatDoNotPass = returnListOfCNVsThatDoNotPass(foundCNVs, bAlleleFreqsNormal[bafsFromThisChr,], bAlleleFreqsTumor[bafsFromThisChr,], 
                                                                  clonalityForChecking, local_purities, local_cn_states, toyBedFile,
                                                                  overdispersionNormal[bafsFromThisChr],
                                                                  overdispersionTumor[bafsFromThisChr],
                                                                  toyLogFoldChange,
                                                                  median(sdsOfProbes) * (sdsOfSomaticSamples[sam_no]),
                                                                  ifelse(sampleInOfftarget, median(sdsOfProbesOff) * (sdsOfSomaticSamplesOff[sam_no_off]), -1)
                                                                  ) 
            if (length(listOfCNVsThatDoNotPass) > 0)
            found_CNVs = found_CNVs[-listOfCNVsThatDoNotPass,,drop=F]
          }
  

          
          if (nrow(found_CNVs) > 0 & !chrom %in% c("chrX", "chrY", "X", "Y")) {
            likeliksFoundCNVsVsPurities <- matrix(0,nrow=nrow(found_CNVs), ncol=length(uniqueLocalPurities))
            for (m in 1:length(uniqueLocalPurities)) {
              localPurityCurrent = uniqueLocalPurities[m]
              for (q in 1:nrow(found_CNVs)) {
                startOfCNV <- found_CNVs[q,2]
                endOfCNV <- found_CNVs[q,3]
                if (endOfCNV - startOfCNV > 3) { 
                  
                  likeliksFoundCNVsVsPurities[q, m] = min(apply(toyMatrixOfLikeliks[(startOfCNV + 1):(endOfCNV - 1),which(local_purities == localPurityCurrent)], 2, sum)
                                                          + penaltyForHigherCN * abs(copy_numbers_for_penalties[which(local_purities == localPurityCurrent)]))
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
                                                  outputFileNameDots, sample_name, toyLogFoldChange, reverseFunctionUsedToTransform, local_cn_states)
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
            cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyLogFoldChange, toyBedFile, output_of_plots, chrom, local_cn_states, local_copy_numbers_used, local_purities, local_cnv_states, 
                                            toySizesOfPointsFromLocalSds, plottingOfPNGs)
            if (found_CNVs[1,1] != -1000) {
              found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
            }
          }
        }
        
      }
      
      
      if (finalIteration == T) {
        if (nrow(found_CNVs_total) > 0){
          makeBarplot(allPotentialPurities, found_CNVs_total)
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
          combinationsOfPurities <- combn(localPurityStates, m)
          #indicesOfPuritiesWithMax = apply(combinationsOfPurities, 2, function(x) {sum(which(x == hundredPercentPurity))})
          #combinationsOfPurities = combinationsOfPurities[,which(indicesOfPuritiesWithMax > 0),drop=F]
          
          bestCombination = 1
          
          for (q in 1:ncol(combinationsOfPurities)) {
            if (length(combinationsOfPurities[,q]) > 1) {
              combinationOfDifferences = combn(uniqueLocalPurities[combinationsOfPurities[,q]], 2)
              allDiffs <- abs(combinationOfDifferences[1,] - combinationOfDifferences[2,])
              #if (min(allDiffs) < 0.0499) {
              #  next
              #}
            }
            minResultForCombination = as.numeric(opt$clonePenalty) * m
            for (r in 1:nrow(likeliksFoundCNVsVsPuritiesGlobal)) {
              minResultForCombination = minResultForCombination + min(
                likeliksFoundCNVsVsPuritiesGlobal[r,combinationsOfPurities[,q]])
            }
            if (minResult > minResultForCombination) {
              bestCombination = q
              resultBestCombination = combinationsOfPurities[,bestCombination]
              minResult = minResultForCombination
            }
          }
        }
        clonalBestPurities = uniqueLocalPurities[resultBestCombination]
        if (length(clonalBestPurities) == 0) {
          clonalBestPurities = c(0, 1)
        }
        
      } else {
        print("No high quality CNVs found in this sample for finding clonality.")
        clonalBestPurities = c(0, 1)
      }
      finalIteration = T
    }
    
    ### STAT TESTS
    CIsOnTarget = matrix(NA, nrow=nrow(found_CNVs_total), ncol=3)
    CIsOnTargetOff = matrix(NA, nrow=nrow(found_CNVs_total), ncol=3)
    BAFsignature = matrix(NA, nrow=nrow(found_CNVs_total), ncol=3)
    overallPvalues = matrix(NA, nrow=nrow(found_CNVs_total), ncol=1)
    if (nrow(found_CNVs_total) > 0){
      for (i in 1:nrow(found_CNVs_total)) {
        defaultCN = 2
        if (found_CNVs_total[i,1] %in% c("chrX", "chrY") & genderOfSamples[germline_sample_no] == "M") {
          defaultCN = 1
        }
        pvalsSeparateTests = c(NA, NA, NA)
        onTargetCoords <- which(bedFile[,1] == found_CNVs_total[i,1] & as.numeric(bedFile[,2]) >= as.numeric(found_CNVs_total[i,2]) & as.numeric(bedFile[,3]) <= as.numeric(found_CNVs_total[i,3]))
        if (length(onTargetCoords) > 1) {
          tumorValue <- median(log2(tumor[onTargetCoords, tumor_sample_no]))
          if (found_CNVs_total[i,1] %in% c("chrX", "chrY")) {
            samplesToUse = which(genderOfSamples == genderOfSamples[germline_sample_no])
          } else {
            samplesToUse = 1:ncol(tmpNormal)
          }
          if (length(samplesToUse) > 2) {
            normalValues <- apply(log2(tmpNormal[onTargetCoords,samplesToUse, drop=F]), 2, median)
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
          offTargetCoords <- which(bedFileOfftarget[,1] == found_CNVs_total[i,1] & as.numeric(bedFileOfftarget[,2]) >= as.numeric(found_CNVs_total[i,2]) & as.numeric(bedFileOfftarget[,3]) <= as.numeric(found_CNVs_total[i,3]))
          if (length(offTargetCoords) > 1) {
            tumorValueOff <- median(log2(tumorOff[offTargetCoords, tumor_sample_no_off]))
            if (found_CNVs_total[i,1] %in% c("chrX", "chrY")) {
              samplesToUseOn = which(genderOfSamples == genderOfSamples[germline_sample_no])
              samplesToUseOff = which(colnames(normalOff) %in% names(samplesToUseOn))
            } else {
              samplesToUseOff = 1:ncol(normalOff)
            }
            if (length(samplesToUseOff) > 2) {
              normalValuesOff <- apply(log2(normalOff[offTargetCoords,samplesToUseOff, drop=F]), 2, median)
              sdOfNormalsOff <- sd(normalValuesOff) * sqrt(matrixWithSdsOff[1, sam_no_off] / median(matrixWithSdsOff[2, ]))
              currentCIOff = c((tumorValueOff), (tumorValueOff + qnorm(0.99) * sdOfNormalsOff), (tumorValueOff + qnorm(0.01) * sdOfNormalsOff))
              currentCIOff = 2 ** (currentCIOff - median(log2(normalOff[offTargetCoords, which(colnames(normalOff) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2])])))
              CIsOnTargetOff[i,] = round(currentCIOff * defaultCN, 2)
            }
            pvalsSeparateTests[2] = 2 * pt( -abs(
              tumorValueOff - median(log2(normalOff[offTargetCoords, which(colnames(normalOff) == strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2])])))
              / sdOfNormalsOff, 
              df=length(samplesToUseOff)) 
          }
        }
        if (!is.null(bAlleleFreqsTumor) & !is.null(bAlleleFreqsNormal)) {
          BAFcoords <- which(bAlleleFreqsNormal[,1] == found_CNVs_total[i,1] & as.numeric(bAlleleFreqsNormal[,2]) >= as.numeric(found_CNVs_total[i,2]) & as.numeric(bAlleleFreqsNormal[,2]) <= as.numeric(found_CNVs_total[i,3]))
          particularAlleleBalance = median(as.numeric(bAlleleFreqsNormal[,5]))
          if (length(BAFcoords) > 0) {
            depthNormal = sum(as.numeric(bAlleleFreqsNormal[BAFcoords,6]))
            normalReversedValues = as.numeric(bAlleleFreqsNormal[BAFcoords,5])
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
      BAFsignature[,3] = p.adjust(as.numeric(BAFsignature[,3]), method="fdr")
      BAFsignature[,3] = format(round(as.numeric(BAFsignature[,3]), 4), scientific = F)
      colnamesForFutureMatrix <- colnames(found_CNVs_total)
      found_CNVs_total = cbind(found_CNVs_total, matrix(CIsOnTarget[,3:2], ncol=2), matrix(CIsOnTargetOff[,3:2], ncol=2), BAFsignature, format(round(overallPvalues,5), scientific = F))
      #colnames(found_CNVs_total) = c(colnamesForFutureMatrix, c("Ontarget_RD", "Ontarget_RD_CI_lower", "Ontarget_RD_CI_upper", "Offtarget_RD", "Offtarget_RD_CI_lower", "Offtarget_RD_CI_upper", "BAF_Normal", "BAF_tumor", "BAF_pval"))
      colnames(found_CNVs_total) = c(colnamesForFutureMatrix, c("Ontarget_RD_CI_lower", "Ontarget_RD_CI_upper", "Offtarget_RD_CI_lower", "Offtarget_RD_CI_upper", "BAF_Normal", "BAF_tumor", "BAF_qval_fdr", "Overall_qvalue"))
    }
    if (length(pvalsForQC > 1)) {
      finalPValue <- 0
    } else {
      finalPValue = 0
    }
    fileToOut <- paste0(folder_name, sample_name, paste0("/CNAs_", sample_name, ".txt"))
    fileConn<-file(fileToOut)
    writeLines(c(paste("##"," QC ", 0, ", gender of sample:", genderOfSamples[germline_sample_no], "clonality by BAF (if != 1):", paste(round(unique(local_purities), digits=3), collapse=";"), collapse = " ")), fileConn)
    close(fileConn)
    if(opt$debug) {
      print(found_CNVs_total)
    }
    found_CNVs_total = found_CNVs_total[order(found_CNVs_total[,1], as.numeric(found_CNVs_total[,2])),]
    write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)	
  }
}

stopCluster(cl)

