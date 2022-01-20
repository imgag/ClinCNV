





cn_states <- c(2,1)
major_alleles = c(1,1)
minor_alleles = c(1,0)
purities <- c(0,0)


purity <- seq(from=0.01 * opt$minimumPurity, to=1, 0.01 * opt$purityStep)

major_allele_cn = 0:30
minor_allele_cn = 0:2

cn_states_tumor_only = c()
actual_major_alleles = c()
actual_minor_alleles = c()

for (pur in purity) {
  for (major in major_allele_cn) {
    for (minor in minor_allele_cn) {
      if (minor < major | (minor == major & (major > 1 | (major == 0 & pur < 0.8)))) {
        if (minor == 0 & major > 4) next
        actual_major_allele = (1 - pur) + pur * major 
        actual_minor_allele = (1 - pur) + pur * minor 
        actual_major_alleles = c(actual_major_alleles, major)
        actual_minor_alleles = c(actual_minor_alleles, minor)
        
        
        cn_states_tumor_only = c(cn_states_tumor_only, actual_major_allele + actual_minor_allele)
        purities <- c(purities, pur)
      }
    }
  }
}

cn_states_tumor_only = c(cn_states, cn_states_tumor_only)
actual_major_alleles = c(major_alleles, actual_major_alleles)
actual_minor_alleles = c(minor_alleles, actual_minor_alleles)


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

folder_name <- paste0(opt$out, "/tumorOnly/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}




covar = F
setwd(opt$out)

setwd(opt$folderWithScript)

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
superRecallThreshold = opt$superRecall


medians_of_coverages <- matrix(nrow=0, ncol=ncol(coverage.normalised))
for (chr in 1:22) {
  which_bed_in_chr = which(bedFileFiltered[,1] == paste0("chr", chr))
  medians_of_coverages <- rbind(medians_of_coverages, apply(coverage.normalised[which_bed_in_chr,], 2, median))
}
for (i in 1:ncol(coverage.normalised)) {
  coverage.normalised[,i] = coverage.normalised[,i] / median(medians_of_coverages[,i])
}
if (frameworkOff == "offtargetGermline") { 
  medians_of_coverages_off <- matrix(nrow=0, ncol=ncol(coverage.normalised.off))
  for (chr in 1:22) {
    which_bed_in_chr = which(bedFileFilteredOfftarget[,1] == paste0("chr", chr))
    medians_of_coverages_off <- rbind(medians_of_coverages_off, apply(coverage.normalised.off[which_bed_in_chr,], 2, median))
  }
  for (i in 1:ncol(coverage.normalised.off)) {
    coverage.normalised.off[,i] = coverage.normalised.off[,i] / median(medians_of_coverages_off[,i])
  }
}


print(paste("Calling started", Sys.time()))
for (sam_no in 1:ncol(coverage.normalised)) {
  outputQCFailed = T
  
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
  
  threshold = opt$scoreS
  minimum_length_of_CNV = opt$lengthS
  price_per_tile = 0.1
  main_initial_state <- 1
  
  
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
  
  
  matrix_of_likeliks_for_FDR <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised), sam_no, localSds, coverage.normalised, sqrt(cn_states_tumor_only / 2))
  
  likeliks_from_baf = NULL
  baf_file = NULL
  # check if the baf file is in the list
  if (!is.null(opt$bafFolder) & file.exists(paste0(opt$bafFolder, "/", sample_name, ".tsv"))) {
    baf_file = read.table(paste0(opt$bafFolder, "/", sample_name, ".tsv"), stringsAsFactors=F, sep="\t")
    baf_file[,5] = as.numeric(baf_file[,5])
    proportion_homozygous_high = c()
    proportion_homozygous_low = c()
    proportion_heterozygous = c(0.48)
    for (chr in unique(baf_file[,1])) {
      if (!chr %in% c("chrX", "chrY")) {
        baf_file_chr = baf_file[which(baf_file[,1] == chr),]
        proportion_homozygous_high <- c(proportion_homozygous_high, 
                                        length(which(baf_file_chr[, 5] > 0.98)) / 
                                          nrow(baf_file_chr))
        
        proportion_homozygous_low <-  c(proportion_homozygous_low, 
                                        length(which(baf_file_chr[, 5] < 0.02)) / 
                                          nrow(baf_file_chr))
        
        proportion_heterozygous <- c(proportion_heterozygous, 
                                     median(baf_file_chr[which(baf_file_chr[, 5] > 0.4 & baf_file_chr[, 5] < 0.6), 5])
        )
      }
    }
    proportion_homozygous_high = median(na.omit(proportion_homozygous_high))
    proportion_homozygous_low = max(median(na.omit(proportion_homozygous_low)), 0.001)
    proportion_heterozygous = median(na.omit(proportion_heterozygous))
    
    
  }
  
  
  if (!is.null(polymorphicRegions)) {
    matrix_of_likeliks_for_FDR[positionsInPolymorphic,] = 0
    matrix_of_likeliks_for_FDR[positionsInPolymorphic,main_initial_state] = 0
  }
  
  fineForTumor = 0.000000001
  matrix_of_likeliks_for_FDR[,which(cn_states_tumor_only %% 1 != 0)] = matrix_of_likeliks_for_FDR[,which(cn_states_tumor_only %% 1 != 0)] + fineForTumor
  matrix_of_likeliks <- matrix_of_likeliks_for_FDR
  
  
  
  if (sam_no_off) {
    matrix_of_likeliks_off = form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised.off), sam_no_off, localSdsOff, coverage.normalised.off, sqrt(cn_states_tumor_only / 2))
    
    globalBed = rbind(bedFileFiltered, bedFileFilteredOfftarget)
    orderOfBed = order(globalBed[,1], as.numeric(globalBed[,2]))
    globalBed = globalBed[orderOfBed,]
    globalMatrixOfLikeliks = rbind(matrix_of_likeliks, matrix_of_likeliks_off)
    globalMatrixOfLikeliks = globalMatrixOfLikeliks[orderOfBed,]
    
    sdsForOutput = c(localSds, localSdsOff)[orderOfBed]
  } else {
    globalBed = bedFileFiltered
    globalMatrixOfLikeliks = matrix_of_likeliks
    
  }
  
  sizesOfPointsFromLocalSds <- 0.1 / localSds 
  
  
  
  maxIteration = 2
  vectorWithNumberOfOutliers <- c()
  vectorOfZScores = rep(0, nrow(bedFileFiltered))
  vectorOfZScores[which(!bedFileFiltered[,1] %in% c("chrX","chrY"))] <- (coverage.normalised[which(!bedFileFiltered[,1] %in% c("chrX","chrY")),sam_no] - 1) / localSds
  if (genderOfSamples[sam_no] == "F") {
    vectorOfZScores[which(bedFileFiltered[,1] %in% c("chrX"))] <- c(vectorOfZScores, (coverage.normalised[which(bedFileFiltered[,1] %in% c("chrX")),sam_no] - 1) / localSds)
  } else {
    vectorOfZScores[which(bedFileFiltered[,1] %in% c("chrX","chrY"))] <- c(vectorOfZScores, (coverage.normalised[which(bedFileFiltered[,1] %in% c("chrX","chrY")),sam_no] - sqrt(1/2)) / localSds)
  }
  
  
  initial_state = main_initial_state
  found_CNVs_total <- matrix(0, nrow=0, ncol=12)
  colnames(found_CNVs_total) <- c("#chr", "start", "end", "CN_change","Purity", "Major allele", "Minor allele", "loglikelihood", "no_of_regions", "length_KB", "potential_AF", "genes")
  
  
  local_cn_states_tumor_only = cn_states_tumor_only
  local_actual_major_alleles = actual_major_alleles
  local_actual_minor_alleles = actual_minor_alleles
  local_purities = purities
  
  matrix_of_likeliks_for_purity = matrix(0, nrow=0, ncol=length(cn_states_tumor_only))
  
  iterations = 0
  iterations = iterations + 1
  
  while (iterations <= 2) {
    if (iterations == 2 & nrow(matrix_of_likeliks_for_purity) > 0) {
      vector_of_allowed_purities <- c(0, 1)
      winning_purities <- vector_of_allowed_purities
      
      matrix_of_likeliks_for_purity_limited = matrix_of_likeliks_for_purity[, which(local_purities %in% vector_of_allowed_purities), drop=F]
      final_sum <- sum(apply(matrix_of_likeliks_for_purity_limited, 1, min))
      
      putential_purities <- unique(local_purities)
      putential_purities = setdiff(putential_purities, vector_of_allowed_purities)
      
      best_purities = list()
      best_purities_scores = list()
      
      for (num_of_clones in 1:4) {
        min_final_sum_cancer = 10**100
        
        penalty_for_additional_clones = as.numeric(opt$clonePenalty) * num_of_clones
        
        combination_of_purities <- combn(putential_purities, num_of_clones)
        winning_combination = combination_of_purities[,1]
        for (i in 1:ncol(combination_of_purities)) {
          vector_of_allowed_purities_cancer <- c(vector_of_allowed_purities, combination_of_purities[,i])
          matrix_of_likeliks_for_purity_limited = matrix_of_likeliks_for_purity[, which(local_purities %in% vector_of_allowed_purities_cancer), drop=F]

          final_sum_cancer <- sum(apply(matrix_of_likeliks_for_purity_limited, 1, min))
          if (final_sum_cancer < min_final_sum_cancer) {
            best_purities[[num_of_clones]] = vector_of_allowed_purities_cancer
            min_final_sum_cancer = final_sum_cancer
            best_purities_scores[[num_of_clones]] = min_final_sum_cancer + penalty_for_additional_clones
          }
        }
      }
      
      for (num_of_clones in 1:4) {
        if (best_purities_scores[[num_of_clones]] < final_sum) {
          winning_purities <- best_purities[[num_of_clones]]
          final_sum = best_purities_scores[[num_of_clones]]
        }
      }
      
      globalMatrixOfLikeliks = globalMatrixOfLikeliks[,which(local_purities %in% winning_purities)]
      
      local_cn_states_tumor_only = local_cn_states_tumor_only[which(local_purities %in% winning_purities)]
      local_actual_major_alleles = local_actual_major_alleles[which(local_purities %in% winning_purities)]
      local_actual_minor_alleles = local_actual_minor_alleles[which(local_purities %in% winning_purities)]
      local_purities = local_purities[which(local_purities %in% winning_purities)]
    }

    

    for (l in 1:length(left_borders)) {
      chrom = names(left_borders)[l]
      if (chrom == "chrX" & genderOfSamples[sam_no] == "M") {
        initial_state <- 2
      } else if (chrom == "chrY" & genderOfSamples[sam_no] == "F") {
        next
      } else if (chrom == "chrY" & genderOfSamples[sam_no] == "M") {
        initial_state <- 2
      } else {
        initial_state <- 1
      }
      
      
      if (!is.null(baf_file) & iterations == 1) {
        baf_from_chr = baf_file[which(baf_file[,1] == chrom & baf_file[,6] > 30),, drop=F]
        if (nrow(baf_from_chr) > 0) {
          
          major_allele = (1 - local_purities) * local_actual_major_alleles[initial_state] + local_purities * local_actual_major_alleles
          minor_allele = (1 - local_purities) * local_actual_major_alleles[initial_state] + local_purities * local_actual_minor_alleles
          
          approximation_baf_expected_value = c(major_allele / (major_allele + minor_allele), minor_allele / (major_allele + minor_allele))
          
          approximation_baf_expected_value = c(0.01, 0.99, approximation_baf_expected_value)
          
          approximation_baf_expected_value = unique(round(approximation_baf_expected_value, 2))
          
          
          likeliks_from_baf = matrix(0, nrow = nrow(baf_from_chr), ncol=length(approximation_baf_expected_value))
          
          for (snv in 1:nrow(baf_from_chr)) {
            corrected_for_alignment = proportion_heterozygous * 2 * approximation_baf_expected_value
            corrected_for_alignment[which.min(corrected_for_alignment)] = 0.01
            corrected_for_alignment[which.max(corrected_for_alignment)] = 0.99
            likeliks_from_baf[snv,] = dbinom(round(baf_from_chr[snv,6] * baf_from_chr[snv,5]), baf_from_chr[snv,6], proportion_heterozygous * 2 * approximation_baf_expected_value )
          }
          likeliks_from_baf[likeliks_from_baf < 10**-50] = 10**-50
          
          actual_likeliks_to_add_chr = matrix(0, nrow = nrow(baf_from_chr), ncol=length(local_cn_states_tumor_only))
          
          hom_loglik_high = likeliks_from_baf[,2]
          hom_loglik_low = likeliks_from_baf[,1]
          
          
          for (column in 1:length(local_cn_states_tumor_only)) {
            higher_baf_index = which(approximation_baf_expected_value == round(major_allele[column] / (major_allele[column] + minor_allele[column]), 2))
            lower_baf_index = which(approximation_baf_expected_value == round(minor_allele[column] / (major_allele[column] + minor_allele[column]), 2))
            
            total_probability = (proportion_homozygous_high * hom_loglik_high +
                                   proportion_homozygous_low * hom_loglik_low + 
                                   (1 - proportion_homozygous_high - proportion_homozygous_low) / 2 * likeliks_from_baf[,higher_baf_index] +
                                   (1 - proportion_homozygous_high - proportion_homozygous_low) / 2 * likeliks_from_baf[,lower_baf_index]
            )
            total_probability = -2 * log(total_probability)
            actual_likeliks_to_add_chr[, column] = total_probability
          }
          
          
          for (snv in 1:nrow(baf_from_chr)) {
            index_to_plus <- which(globalBed[,1] == chrom & globalBed[,2]  <= baf_from_chr[snv,2] & globalBed[,3]  >= baf_from_chr[snv,3])
            if (length(index_to_plus) > 0) {
              globalMatrixOfLikeliks[index_to_plus,] = globalMatrixOfLikeliks[index_to_plus,] + actual_likeliks_to_add_chr[snv,]
            }
          }
        }
        
      }
      
      
      
      
      
      start = left_borders[[l]]
      end = right_borders[[l]]
      for (k in 1:2) {
        local_cn_states = local_cn_states_tumor_only
        
        
        output_of_plots <-  paste0(folder_name, sample_name)
        which_to_allow <- "NA"
        which_to_allow_ontarget <- "NA"
        startOfChrom = 0
        if (chrom == "chrX" & !is.na(startX)) {
          startOfChrom = startX
        }
        if (chrom == "chrXP") {
          chrom = "chrX"
        }
        if (k == 1) {
          which_to_allow = which(globalBed[,1] == chrom & as.numeric(globalBed[,3]) <= as.numeric(left_borders[[l]]) &  as.numeric(globalBed[,2]) >= startOfChrom)
          which_to_allow_ontarget = which(bedFileFiltered[,1] == chrom & as.numeric(bedFileFiltered[,3]) <= as.numeric(left_borders[[l]]) &  as.numeric(bedFileFiltered[,2]) >= startOfChrom)
        } else {
          which_to_allow = which(globalBed[,1] == chrom & as.numeric(globalBed[,2]) >= as.numeric(right_borders[[l]]) & as.numeric(globalBed[,3]) <= ends_of_chroms[[l]])
          which_to_allow_ontarget = which(bedFileFiltered[,1] == chrom & as.numeric(bedFileFiltered[,2]) >= as.numeric(right_borders[[l]]) & as.numeric(bedFileFiltered[,3]) <= ends_of_chroms[[l]] )
        }
        
        
        toyMatrixOfLikeliks = globalMatrixOfLikeliks[which_to_allow,]
        toybedFileFiltered = globalBed[which_to_allow,]
        
        copy_numbers_for_penalties = 3 - (local_actual_major_alleles + local_actual_minor_alleles)
        copy_numbers_for_penalties[which(copy_numbers_for_penalties > 0)] = 0
        toyMatrixOfLikeliks = sweep(toyMatrixOfLikeliks, 2, abs(copy_numbers_for_penalties) * opt$pnealtyHigherCopyOneSegment, FUN="+")
        
        penaltyForHigherCN = opt$pnealtyHigherCopy
        copy_numbers_for_penalties = 3 - (local_actual_major_alleles + local_actual_minor_alleles)
        copy_numbers_for_penalties[which(copy_numbers_for_penalties > 0)] = 0
        penalties = penaltyForHigherCN * abs(copy_numbers_for_penalties)
        
        found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, initial_state, c(), penalties))
        
        
        toySizesOfPointsFromLocalSds = sizesOfPointsFromLocalSds[which_to_allow]
        toyCoverageGermline = coverage.normalised[which_to_allow_ontarget,sam_no]
        toyCoverageGermlineCohort = coverage.normalised[which_to_allow_ontarget,]
        
        if (nrow(found_CNVs) > 0 & iterations == 1) {
          for (cnv in 1:nrow(found_CNVs)) {
            if (found_CNVs[,1] < -2 * threshold) {
              likeliks_to_decide = apply(toyMatrixOfLikeliks[found_CNVs[cnv, 2]:found_CNVs[cnv, 3],], 2, sum)
              
              matrix_of_likeliks_for_purity = rbind(matrix_of_likeliks_for_purity, likeliks_to_decide)
            }
          }
        }
        
        
        
        
        if (nrow(found_CNVs) > 0 & iterations == 2) {
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
        
        if (!chrom %in% c("chrX", "chrY") & iterations == 2) {
          vectorOfZScoresLocaL = vectorOfZScores[which_to_allow]
          if (length(vectorOfZScoresLocaL) > 10)
            vectorWithNumberOfOutliers = c(vectorWithNumberOfOutliers, 
                                           length(which(vectorOfZScoresLocaL > qnorm(0.975) | vectorOfZScoresLocaL < qnorm(0.025))) / length(vectorOfZScoresLocaL))
        }
        
        ### IGV PLOTTING
        if (opt$visulizationIGV & iterations == 2) {
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
                                              outputFileNameDots, sample_name, toyCoverageGermline, reverseFunctionUsedToTransform, local_cn_states, sdsForOutput[which_to_allow], outputQCFailed)
          outputQCFailed = F
          if(opt$debug) {
            print("END OF IGV PLOTTING")
          }
        }
        ### END OF IGV PLOTTING
        
        
        
        if (nrow(found_CNVs) > 0 & iterations == 2) {
          # UNCOMMENT FOR PLOTTING!!!
          
          cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyCoverageGermline, toybedFileFiltered, output_of_plots, chrom, local_cn_states, 
                                          toySizesOfPointsFromLocalSds,alleleFrequency, plottingOfPNGs)
          
          vec_of_purities <- local_purities[found_CNVs[,4]]
          vec_of_major_alleles <- local_actual_major_alleles[found_CNVs[,4]]
          vec_of_minor_alleles <- local_actual_minor_alleles[found_CNVs[,4]]

          
          cnvsToWriteOut_tmp = matrix(nrow=nrow(cnvsToWriteOut), ncol=ncol(found_CNVs_total))
          cnvsToWriteOut_tmp[,1:4] = cnvsToWriteOut[,1:4]
          
          cnvsToWriteOut_tmp[,8:ncol(found_CNVs_total)] = cnvsToWriteOut[,5:ncol(cnvsToWriteOut)]
          cnvsToWriteOut_tmp[,5] = vec_of_purities
          cnvsToWriteOut_tmp[,6] = vec_of_major_alleles
          cnvsToWriteOut_tmp[,7] = vec_of_minor_alleles
          
          if (found_CNVs[1,1] != -1000) {
            found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut_tmp)
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
    iterations = iterations + 1
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
              genderOfThisSample = genderOfSamples[which(names(genderOfSamples) == sample_name)]
              idsOfSameSexSamples = which(colnames(coverage.normalised.off) %in% names(genderOfSamples)[which(genderOfSamples == genderOfThisSample)])
              valueOfOthers = apply(coverage.normalised.off[coordsInBedOff,idsOfSameSexSamples,drop=F],2,median)
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
  found_CNVs_total[,10] = (format(as.numeric(found_CNVs_total[,10]), nsmall=3))
  found_CNVs_total[,11] = (format(as.numeric(found_CNVs_total[,11]), nsmall=3))
  found_CNVs_total[which(found_CNVs_total[,1] == "chrXP"),1] = "chrX"
  found_CNVs_total[which(is.na(found_CNVs_total[,12]) | found_CNVs_total[,12] == "na"),12] = "" 
  write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)
}


