library(RColorBrewer)


prepareDataAndCall <- function(bedFileForCluster, tmpNormal, tumor, genderOfSamples, bedFileForClusterOff=NULL, tmpNormalOff=NULL, tumorOff=NULL) {
  no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
  cl<-makeCluster(no_cores, type="FORK")
  registerDoParallel(cl)
  
  ### PROCESSING OF SOMATIC VARIANTS
  setwd(opt$folderWithScript)
  source(paste0(opt$folderWithScript, "/somatic/helpersSomatic.R"),local=TRUE)
  source(paste0(opt$folderWithScript, "/somatic/helpersBalleleFreq.R"),local=T)
  source(paste0(opt$folderWithScript, "/somatic/somaticCallingProcedure.R"),local=T)
  
  
  
  
  print(paste("Work on data preparation for somatic samples started (log-fold change matrices plus parameters estimation)", Sys.time()))
  listOfTmpNormalAndTmpTumor = normalizeToCommonMedian(tmpNormal, tumor, bedFileForCluster, genderOfSamples)
  tmpNormal = listOfTmpNormalAndTmpTumor[[1]] 
  tmpTumor = listOfTmpNormalAndTmpTumor[[2]] 
  rm(listOfTmpNormalAndTmpTumor)
  
  listOfValue <- formilngLogFoldChange(pairs, tmpNormal, tmpTumor, bedFileForCluster, genderOfSamples)
  matrixOfLogFold <- listOfValue[[1]]
  gendersInOntargetMatrix = listOfValue[[2]]
  
  
  
  matrixWithSdsList <- findSDsOfSamples(pairs, tmpNormal, tmpTumor, bedFileForCluster, bordersOfChroms, genderOfSamples)
  matrixWithSds = matrixWithSdsList[[1]]
  sdsOfSomaticSamples <- matrixWithSds[4,]
  sdsOfProbes <- matrixWithSdsList[[2]]
  
  # QC control
  sampleLevelOfNoise = apply(matrixOfLogFold[which(!bedFileForCluster[,1] %in% c("chrX","chrY")),], 2, Sn)
  samplesNotPassedQC = colnames(matrixOfLogFold)[which(sampleLevelOfNoise > 1.0)]
  
  if (length(samplesNotPassedQC) > 0) {
    print("Samples don't pass our QC - too noisy!")
    print(samplesNotPassedQC)
    if (!is.null(opt$normalSample) & !is.null(opt$tumorSample)) {
      if (paste(opt$tumorSample, opt$normalSample, sep="-") %in% samplesNotPassedQC) {
        print(paste("But you explicitly specifed the sample name as", paste(opt$tumorSample, opt$normalSample, sep="-"), "so we still proceed with its analysis"))
        samplesNotPassedQC = setdiff(samplesNotPassedQC, paste(opt$tumorSample, opt$normalSample, sep="-"))
      }
    }
    if (length(samplesNotPassedQC) > 0) {
      sdsOfSomaticSamples = sdsOfSomaticSamples[-which(colnames(matrixOfLogFold) %in% samplesNotPassedQC)]
      matrixWithSds = matrixWithSds[,-which(colnames(matrixOfLogFold) %in% samplesNotPassedQC)]
      gendersInOntargetMatrix = gendersInOntargetMatrix[-which(colnames(matrixOfLogFold) %in% samplesNotPassedQC)]
      matrixOfLogFold = matrixOfLogFold[,-which(colnames(matrixOfLogFold) %in% samplesNotPassedQC)]
    }
  }
  
  ## QC FOR PROBES
  listOfProbesQC = probeLevelQC(matrixOfLogFold, sdsOfProbes, sdsOfSomaticSamples, gendersInOntargetMatrix, bedFileForCluster)
  probesToRemove <- listOfProbesQC[[1]]
  sdsOfProbes = listOfProbesQC[[2]]
  if (length(probesToRemove) > 0) {
    sdsOfProbes = sdsOfProbes[-probesToRemove]
    matrixOfLogFold = matrixOfLogFold[-probesToRemove,]
    bedFileForCluster = bedFileForCluster[-probesToRemove,]
    tmpNormal = tmpNormal[-probesToRemove,]
    tmpTumor = tmpTumor[-probesToRemove,]
  }
  
  if (frameworkOff == "offtarget") {
    listOfTmpNormalAndTmpTumor = normalizeToCommonMedian(tmpNormalOff, tumorOff, bedFileForClusterOff, genderOfSamples)
    tmpNormalOff = listOfTmpNormalAndTmpTumor[[1]]
    tmpTumorOff = listOfTmpNormalAndTmpTumor[[2]]
    listOfValueOff <- formilngLogFoldChange(pairs, tmpNormalOff[,which(colnames(tmpNormalOff) %in% colnames(tmpNormal))], tmpTumorOff, bedFileForClusterOff, genderOfSamples)
    matrixOfLogFoldOff =  listOfValueOff[[1]]
    gendersInOfftargetMatrix = listOfValueOff[[2]]
    
    if (length(samplesNotPassedQC) > 0) {
      if (length(which(colnames(matrixOfLogFoldOff) %in% samplesNotPassedQC)) > 0) {
        gendersInOfftargetMatrix = gendersInOfftargetMatrix[-which(colnames(matrixOfLogFoldOff) %in% samplesNotPassedQC)]
        matrixOfLogFoldOff = matrixOfLogFoldOff[,-which(colnames(matrixOfLogFoldOff) %in% samplesNotPassedQC)]
      }
    }
    
    
    
    
    bordersOfChroms <- getBordersOfChromosomes(bedFileForClusterOff)
    matrixWithSdsOffList <- findSDsOfSamples(pairs, tmpNormalOff[,which(colnames(tmpNormalOff) %in% colnames(tmpNormal))], tmpTumorOff, bedFileForClusterOff, bordersOfChroms, genderOfSamples)
    matrixWithSdsOff = matrixWithSdsOffList[[1]]
    sdsOfSomaticSamplesOff <- matrixWithSdsOff[4,]
    sdsOfProbesOff <- matrixWithSdsOffList[[2]]
    
    listOfProbesQC <- probeLevelQC(matrixOfLogFoldOff, sdsOfProbesOff, sdsOfSomaticSamplesOff, gendersInOfftargetMatrix, bedFileForClusterOff)
    probesToRemove <- listOfProbesQC[[1]]
    sdsOfProbesOff = listOfProbesQC[[2]]
    if (length(probesToRemove) > 0) {
      sdsOfProbesOff = sdsOfProbesOff[-probesToRemove]
      matrixOfLogFoldOff = matrixOfLogFoldOff[-probesToRemove,]
      bedFileForClusterOff = bedFileForClusterOff[-probesToRemove,]
      tmpNormalOff = tmpNormalOff[-probesToRemove,]
      tmpTumorOff = tmpTumorOff[-probesToRemove,]
    }
    
  }
  
  
  
  
  
  
  cn_states <- c()
  copy_numbers = 0:15
  purity <- seq(from=opt$minimumPurity, to=100.1, by=opt$purityStep) / 100
  purities <- c()
  copy_numbers_used_major = c()
  copy_numbers_used_minor = c()
  BAF_number_of_reads_minor = c()
  BAF_number_of_reads_major = c()
  ### DESCRIPTION OF STATES
  # CNV - copy number change, 1 allele changed
  # CNVboth - duplication when both alleles changed
  # LOH - Loss of Heterozygosity
  # normal - nothing changed comparing to normal genome
  # CNVcomplex - not single allelic CNV
  
  for (pur in purity) {
    for (cn in copy_numbers) {
      for (cn1 in copy_numbers) {
        if (cn1 <= 4 & cn1 <= cn & !(cn1 == 1 & cn == 1) & 
            (cn1 + cn) <= max(copy_numbers) & ((1-pur) + pur * cn) > 0 ) {
          if (cn1 == cn & pur < opt$clonalityForChecking) next
          #if (cn1 == 0 & cn > 8) next
          #if (cn1 == cn & cn1 + cn > 8) next
          cn_state = (1 - pur) * 2 + pur * cn + pur * cn1
          cn_states <- c(cn_states, cn_state)
          purities <- c(purities, pur)
          copy_numbers_used_major <- c(copy_numbers_used_major, cn)
          copy_numbers_used_minor <- c(copy_numbers_used_minor, cn1)
          BAF_number_of_reads_minor <- c(BAF_number_of_reads_minor, (1-pur) + pur * cn1)
          BAF_number_of_reads_major <- c(BAF_number_of_reads_major, (1-pur) + pur * cn)
        }
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
  
  datasetOfPuritiesCopies <- cbind(cn_states, copy_numbers_used_major, copy_numbers_used_minor, BAF_number_of_reads_major, BAF_number_of_reads_minor, purities)
  colnames(datasetOfPuritiesCopies) <- c("state", "major_cn", "minor_cn", "major_BAF", "minor_BAF", "puritiy")
  
  final_order <- order(cn_states)
  datasetOfPuritiesCopies = datasetOfPuritiesCopies[final_order,]
  datasetOfPuritiesCopies = rbind(c(2,1,1,1,1,0), c(1,1,0,1,0,0), datasetOfPuritiesCopies)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  folder_name <- paste0(opt$out, "/somatic/")
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  allPotentialPurities <- unique(purities)
  penaltyForHigherCN = 5
  penaltyForHigherCNoneTile = 0.05
  clonalityForChecking = opt$clonalityForChecking
  print(paste("Work on actual calling started.", Sys.time()))
  
  if (frameworkDataTypes == "covdepthBAF") {
    normalNames = sapply(1:length(allowedChromsBaf), function(i) {strsplit(names(allowedChromsBaf)[i], split="-")[[1]][2]})
  } else {
    normalNames = c()
  }
  
  
  somaticCalling(matrixOfLogFold)
  stopCluster(cl)
}




