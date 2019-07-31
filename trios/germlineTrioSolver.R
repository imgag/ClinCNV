

cn_states <- 0:6


vect_of_norm_likeliks <- fast_dt_list(as.numeric(opt$degreesOfFreedomStudent))



startCoordOfNonInterruptedSegment = 1
shiftsOfCoverage <- c()
colours = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414,30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414)]

overallResult <- matrix(0, nrow=0, ncol=7)
folder_name <- paste0(opt$out, "/trios/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

vectors_of_cn_states = matrix(c(2,2,2), nrow=1,ncol=3)
# ADD DE NOVO CNVs
for (i in 1:length(cn_states)) {
  vectors_of_cn_states = rbind(vectors_of_cn_states, matrix(c(cn_states[i],2,2), nrow=1,ncol=3))
  vectors_of_cn_states = rbind(vectors_of_cn_states, matrix(c(cn_states[i],2,1), nrow=1,ncol=3)) # X chromosome
  vectors_of_cn_states = rbind(vectors_of_cn_states, matrix(c(cn_states[i],0,1), nrow=1,ncol=3)) # Y chromosome
}
# ADD COMBINATIONS OF PARENTS CNVs
for (i in 1:length(cn_states)) {
  for (j in 1:length(cn_states)) {
    startForFathers = ifelse(cn_states[i] >= 2, 1, 0)
    endForFathers = ifelse(cn_states[i] >= 2, cn_states[i] - 1, cn_states[i])
    father_alleles_1 = startForFathers:endForFathers
    father_alleles_2 = cn_states[i] - father_alleles_1
    startForMothers = ifelse(cn_states[j] >= 2, 1, 0)
    endForMothers = ifelse(cn_states[j] >= 2, cn_states[j] - 1, cn_states[j])
    
    mother_alleles_1 = startForMothers:endForMothers
    mother_alleles_2 = cn_states[j] - mother_alleles_1
    for (k in 1:length(father_alleles_1)) {
      for (l in 1:length(mother_alleles_1)) { 
        vectors_of_cn_states = rbind(vectors_of_cn_states, matrix(c(father_alleles_1[k] + mother_alleles_1[l], cn_states[i],cn_states[j]), nrow=1,ncol=3))
        vectors_of_cn_states = rbind(vectors_of_cn_states, matrix(c(father_alleles_1[k] + mother_alleles_2[l], cn_states[i],cn_states[j]), nrow=1,ncol=3))
        vectors_of_cn_states = rbind(vectors_of_cn_states, matrix(c(father_alleles_2[k] + mother_alleles_1[l], cn_states[i],cn_states[j]), nrow=1,ncol=3))
        vectors_of_cn_states = rbind(vectors_of_cn_states, matrix(c(father_alleles_2[k] + mother_alleles_2[l], cn_states[i],cn_states[j]), nrow=1,ncol=3))
      }
    }
  }
}
vectors_of_cn_states = unique(vectors_of_cn_states)
#vectors_of_cn_states = vectors_of_cn_states[which(vectors_of_cn_states[,1] != 2 | (vectors_of_cn_states[,1] == 2 & vectors_of_cn_states[,2] == 2 & vectors_of_cn_states[,3] == 2)),]


setwd(opt$folderWithScript)
source("./trios/helpersTrio.R")

for (trioRow in 1:nrow(trios)) {
  child_number <- which(colnames(coverage.normalised) == trios[trioRow, 1])
  mother_number  <- which(colnames(coverage.normalised) == trios[trioRow, 2])
  father_number  <- which(colnames(coverage.normalised) == trios[trioRow, 3])

  
  sample_name = paste(trios[trioRow,], collapse="-")
  
  if (length(child_number) == 0 | length(father_number) == 0 | length(mother_number) == 0) {
    print("Trio with names")
    print(trios[trioRow,])
    print("does not exist, we try next")
    next()
  } else {
    print("We analyse trio")
    print(trios[trioRow,])
    child_name = colnames(coverage.normalised)[child_number]
    mother_name = colnames(coverage.normalised)[mother_number]
    father_name = colnames(coverage.normalised)[father_number]
  }
  
  threshold = opt$scoreG
  minimum_length_of_CNV = opt$lengthG
  price_per_tile = 1
  initial_state <- 1
  
  
  localSdsChild = sdsOfProbes * sdsOfGermlineSamples[child_number] 
  localSdsMother = sdsOfProbes * sdsOfGermlineSamples[mother_number] 
  localSdsFather = sdsOfProbes * sdsOfGermlineSamples[father_number] 
  
  if (!dir.exists(paste0(folder_name, sample_name))) {
    dir.create(paste0(folder_name, sample_name))
  } else {
    if (!is.null(opt$reanalyseCohort)) {
      next
    }
  }
  setwd(paste0(folder_name, sample_name))
  
  
  dict_to_output = c()
  
  cn_states_child = sort(unique(vectors_of_cn_states[,1]))
  cn_states_mother = sort(unique(vectors_of_cn_states[,2]))
  cn_states_father = sort(unique(vectors_of_cn_states[,3]))
  
  matrix_of_likeliks_child <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised), child_number, localSdsChild, coverage.normalised, sqrt(cn_states_child / 2))
  matrix_of_likeliks_mother <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised), mother_number, localSdsMother, coverage.normalised, sqrt(cn_states_mother / 2))
  matrix_of_likeliks_father <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised), father_number, localSdsFather, coverage.normalised, sqrt(cn_states_father / 2))
  
  globalBed = bedFile
  matrix_of_likeliks_child_full = matrix_of_likeliks_child
  matrix_of_likeliks_mother_full = matrix_of_likeliks_mother
  matrix_of_likeliks_father_full = matrix_of_likeliks_father
  coverageChildFull = coverage.normalised[,child_number]
  coverageMotherFull = coverage.normalised[,mother_number]
  coverageFatherFull = coverage.normalised[,father_number]
  if (frameworkOff == "offtargetGermline") {
    child_number_off <- which(colnames(coverage.normalised.off) == trios[trioRow, 1])
    mother_number_off  <- which(colnames(coverage.normalised.off) == trios[trioRow, 2])
    father_number_off  <- which(colnames(coverage.normalised.off) == trios[trioRow, 3])
    if (length(child_number_off) == 1 | length(father_number_off) == 1 | length(mother_number_off) == 1) {
      localSdsChildOff = sdsOfProbesOff * sdsOfGermlineSamplesOff[child_number] 
      localSdsMotherOff = sdsOfProbesOff * sdsOfGermlineSamplesOff[mother_number] 
      localSdsFatherOff = sdsOfProbesOff * sdsOfGermlineSamplesOff[father_number] 
      matrix_of_likeliks_child_off <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised.off), child_number, localSdsChild, coverage.normalised.off, sqrt(cn_states_child / 2))
      matrix_of_likeliks_mother_off <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised.off), mother_number, localSdsMother, coverage.normalised.off, sqrt(cn_states_mother / 2))
      matrix_of_likeliks_father_off <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised.off), father_number, localSdsFather, coverage.normalised.off, sqrt(cn_states_father / 2))
      
      matrix_of_likeliks_child_full = rbind(matrix_of_likeliks_child, matrix_of_likeliks_child_off)
      matrix_of_likeliks_mother_full = rbind(matrix_of_likeliks_mother, matrix_of_likeliks_mother_off)
      matrix_of_likeliks_father_full = rbind(matrix_of_likeliks_father, matrix_of_likeliks_father_off)
      
      globalBed = rbind(bedFile, bedFileOfftarget)
      
      orderOfGlobalBed = order(globalBed[,1], as.numeric(globalBed[,2]))
      matrix_of_likeliks_child_full = matrix_of_likeliks_child_full[orderOfGlobalBed,]
      matrix_of_likeliks_mother_full = matrix_of_likeliks_mother_full[orderOfGlobalBed,]
      matrix_of_likeliks_father_full = matrix_of_likeliks_father_full[orderOfGlobalBed,]
      globalBed = globalBed[orderOfGlobalBed,]
      
      coverageChildFull = c(coverageChildFull, coverage.normalised.off[,child_number_off] )
      coverageMotherFull = c(coverageMotherFull, coverage.normalised.off[,mother_number_off] )
      coverageFatherFull = c(coverageFatherFull, coverage.normalised.off[,father_number_off] )
      coverageChildFull = coverageChildFull[orderOfGlobalBed]
      coverageMotherFull = coverageMotherFull[orderOfGlobalBed]
      coverageFatherFull = coverageFatherFull[orderOfGlobalBed]
    }
  } 

  
  matrix_of_likeliks = matrix(0, nrow=nrow(matrix_of_likeliks_child_full), ncol=nrow(vectors_of_cn_states))
  for (i in 1:nrow(vectors_of_cn_states)) {
    child_state = which(cn_states_child == vectors_of_cn_states[i,1])
    mother_state = which(cn_states_mother == vectors_of_cn_states[i,2])
    father_state = which(cn_states_father == vectors_of_cn_states[i,3])
    
    matrix_of_likeliks[,i] = matrix_of_likeliks_child_full[,child_state] + matrix_of_likeliks_mother_full[,mother_state] + matrix_of_likeliks_father_full[,father_state]
  }
  
  matrix_of_likeliks_read_depth_only = matrix_of_likeliks
  
  sizesOfPointsFromLocalSds <- 0.1 / localSdsChild 
  
  
  
  numberOfCNVsIsSufficientlySmall = F
  iterations = 0
  maxIteration = opt$maxNumIter
  while (!numberOfCNVsIsSufficientlySmall & iterations < maxIteration) {
    found_CNVs_total <- matrix(0, nrow=0, ncol=14)
    colnames(found_CNVs_total) <- c("#chr", "start", "end", "kid_CN_change", "mother_CN_change", "father_CN_change", "priority", "loglikelihood", "no_of_regions", "length_in_KB", "loglik_kid", "loglik_mother", "loglik_father", "genes")
    
    iterations = iterations + 1
    for (l in 1:length(left_borders)) {
      chrom = names(left_borders)[l]
      start = left_borders[[l]]
      end = right_borders[[l]]
      for (k in 1:2) {
        if (nrow(found_CNVs_total) > opt$maxNumGermCNVs * 3) {
          break
        }
        local_vectors_of_cn_states = vectors_of_cn_states
        if (chrom %in% c("chrX", "X")) {
          if (genderOfSamples[child_number] == "F")
            local_vectors_of_cn_states[1,] = c(2,2,1)
          else 
            local_vectors_of_cn_states[1,] = c(1,2,1)
        }
        if (chrom %in% c("chrY", "Y")) {
          if (genderOfSamples[child_number] == "F")
          local_vectors_of_cn_states[1,] = c(0,0,1)
          else
            local_vectors_of_cn_states[1,] = c(1,0,1)
        }
        output_of_plots <-  paste0(folder_name, sample_name)
        which_to_allow <- "NA"
        if (k == 1) {
          which_to_allow = which(globalBed[,1] == chrom & globalBed[,2] <= as.numeric(start) )
        } else {
          which_to_allow = which(globalBed[,1] == chrom & globalBed[,2] >= as.numeric(end) )
        }
        toyMatrixOfLikeliks = matrix_of_likeliks[which_to_allow,]
        toyBedFile = globalBed[which_to_allow,]
        found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, initial_state))
        
        

        
        ### IGV PLOTTING
        
        toyLogFoldChangeChild = coverageChildFull[which_to_allow]
        toySizesOfPointsFromLocalSdsChild = sizesOfPointsFromLocalSds[which_to_allow]
        outputFileNameCNVsChild <- paste0(folder_name, sample_name, "/", child_name, "_cnvs.seg")
        outputFileNameDotsChild <- paste0(folder_name, sample_name, "/", child_name, "_cov.seg")
        reverseFunctionUsedToTransform = function(x, chrom) {return((2 * x ** 2))}
        outputSegmentsAndDotsFromListOfCNVs(toyBedFile, found_CNVs, start, end, outputFileNameCNVsChild, 
                                            outputFileNameDotsChild, paste0(child_name, "_kid"), toyLogFoldChangeChild, reverseFunctionUsedToTransform, local_vectors_of_cn_states[,1])
        
        toyLogFoldChangeMother = coverageMotherFull[which_to_allow]
        toySizesOfPointsFromLocalSdsMother = sizesOfPointsFromLocalSds[which_to_allow]
        outputFileNameCNVsMother <- paste0(folder_name, sample_name, "/", mother_name, "_cnvs.seg")
        outputFileNameDotsMother <- paste0(folder_name, sample_name, "/", mother_name, "_cov.seg")
        reverseFunctionUsedToTransform = function(x, chrom) {return((2 * x ** 2))}
        outputSegmentsAndDotsFromListOfCNVs(toyBedFile, found_CNVs, start, end, outputFileNameCNVsMother, 
                                            outputFileNameDotsMother, paste0(mother_name, "_mother"), toyLogFoldChangeMother, reverseFunctionUsedToTransform, local_vectors_of_cn_states[,2])
        
        toyLogFoldChangeFather = coverageFatherFull[which_to_allow]
        toySizesOfPointsFromLocalSdsFather = sizesOfPointsFromLocalSds[which_to_allow]
        outputFileNameCNVsFather <- paste0(folder_name, sample_name, "/", father_name, "_cnvs.seg")
        outputFileNameDotsFather <- paste0(folder_name, sample_name, "/", father_name, "_cov.seg")
        reverseFunctionUsedToTransform = function(x, chrom) {return((2 * x ** 2))}
        outputSegmentsAndDotsFromListOfCNVs(toyBedFile, found_CNVs, start, end, outputFileNameCNVsFather, 
                                            outputFileNameDotsFather, paste0(father_name, "_father"), toyLogFoldChangeFather, reverseFunctionUsedToTransform, local_vectors_of_cn_states[,3])
        
        if(opt$debug) {
          print("END OF IGV PLOTTING")
        }
        ### END OF IGV PLOTTING
        
        
        
        if (nrow(found_CNVs) > 0) {
          # UNCOMMENT FOR PLOTTING!!!
          
          matrixOfScoresSonMomFather <- matrix(0, ncol=3, nrow=0)
          for (s in 1:nrow(found_CNVs)) {
            statesOfCNVsInTrio = local_vectors_of_cn_states[found_CNVs[s,4],]
            likeliksKid = sum(matrix_of_likeliks_child_full[which_to_allow,][found_CNVs[s,2]:found_CNVs[s,3],statesOfCNVsInTrio[1] + 1]) - sum(matrix_of_likeliks_child_full[which_to_allow,][found_CNVs[s,2]:found_CNVs[s,3],local_vectors_of_cn_states[1,1] + 1])
            likeliksMadre = sum(matrix_of_likeliks_mother_full[which_to_allow,][found_CNVs[s,2]:found_CNVs[s,3],statesOfCNVsInTrio[2] + 1]) - sum(matrix_of_likeliks_mother_full[which_to_allow,][found_CNVs[s,2]:found_CNVs[s,3],local_vectors_of_cn_states[1,2] + 1])
            likeliksPadre = sum(matrix_of_likeliks_father_full[which_to_allow,][found_CNVs[s,2]:found_CNVs[s,3],statesOfCNVsInTrio[3] + 1]) - sum(matrix_of_likeliks_father_full[which_to_allow,][found_CNVs[s,2]:found_CNVs[s,3],local_vectors_of_cn_states[1,3] + 1])
            matrixOfScoresSonMomFather = rbind(matrixOfScoresSonMomFather, matrix(c(-1 * round(likeliksKid, 2), round(-1 * likeliksMadre, 2),  round(-1 * likeliksPadre, 2)), nrow=1))
          }
          
          cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyLogFoldChange, toyBedFile, output_of_plots, chrom, local_vectors_of_cn_states, 
                                          toySizesOfPointsFromLocalSds, matrixOfScoresSonMomFather, plottingOfPNGs)
          if (found_CNVs[1,1] != -1000) {
            found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
            if (nrow(found_CNVs_total) > 3 * opt$maxNumGermCNVs) {
              break
            }
          }
          for (i in 1:nrow(found_CNVs)) {
            
            CNVnamesInside <- unlist(unique(toyBedFile[found_CNVs[i,2]:found_CNVs[i,3],4]))
            if(opt$debug) {
              print(CNVnamesInside)
            }
            
            CNVentry = matrix(c(sample_name, chrom, toyBedFile[found_CNVs[i,2],2], toyBedFile[found_CNVs[i,3],3], 
                                paste(CNVnamesInside, collapse=", "),
                                found_CNVs[i,4] - 1, 
                                found_CNVs[i,5]),
                              nrow=1)
            if(opt$debug) {
              print(CNVentry)
            }
            overallResult = rbind(overallResult, CNVentry)
          }
        }
      }
      if (nrow(found_CNVs_total) > 3 * opt$maxNumGermCNVs & iterations != maxIteration) {
        break
      }
    }
    if (nrow(found_CNVs_total) < opt$maxNumGermCNVs * 3) {
      numberOfCNVsIsSufficientlySmall = T
      iterations = 4
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
  finalPValue = 1.0
  fileToOut <- paste0(folder_name, sample_name, paste0("/", sample_name, "_CNVs.txt"))
  fileConn<-file(fileToOut)
  writeLines(c(paste("##"," QC ", finalPValue, collapse = " ")), fileConn)
  close(fileConn)
  found_CNVs_total = found_CNVs_total[order(as.numeric(found_CNVs_total[,7]), as.numeric(found_CNVs_total[,8]), decreasing = T), ]
  write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)
}

