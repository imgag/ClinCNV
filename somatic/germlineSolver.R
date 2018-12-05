autosomes <- which(!bedFile[,1] %in% c("chrX", "chrY", "X", "Y"))
sdsOfGermlineSamples <- apply(coverage.normalised[autosomes,], 2, determineSDsOfGermlineSample)





#locationsShiftedLogFoldChanges <- sweep(matrixOfLogFold, 1, locations)

listOfVarianceAndMultiplicator <- esimtateVarianceFromSampleNoise(sdsOfGermlineSamples, 100000)
esimtatedVarianceFromSampleNoise <- listOfVarianceAndMultiplicator[[2]]
multiplicator <- listOfVarianceAndMultiplicator[[1]]

vect_of_t_likeliks <- fast_dt_list(ncol(coverage.normalised) - 1)
vect_of_norm_likeliks <- fast_dnorm_list()


cn_states <- 0:20








startCoordOfNonInterruptedSegment = 1
shiftsOfCoverage <- c()
colours = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414,30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414)]

overallResult <- matrix(0, nrow=0, ncol=7)
folder_name <- paste0(opt$out, "/normal/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}


for (sam_no in 1:ncol(coverage.normalised)) {
  sample_name <- colnames(coverage.normalised)[sam_no]
  
  if (!is.null(opt$normalSample)) {
    if (!sample_name == opt$normalSample) {
      next
    }
  }
  
  threshold = opt$scoreG
  minimum_length_of_CNV = opt$lengthG
  price_per_tile = 1
  initial_state <- 3
  
  
  localSds = sdsOfProbes * esimtatedVarianceFromSampleNoise[sam_no] * multiplicator
  
  
  
  
  if(opt$debug) {
    print(sam_no)
  }
  if(opt$debug) {
    print(sample_name)
  }
  if (!dir.exists(paste0(folder_name, sample_name))) {
    dir.create(paste0(folder_name, sample_name))
  }
  setwd(paste0(folder_name, sample_name))
  
  
  dict_to_output = c()
  
  
  
  
  matrix_of_likeliks <- form_matrix_of_likeliks_one_sample(1, ncol(coverage.normalised), sam_no, localSds, coverage.normalised, sqrt(cn_states / 2))
  matrix_of_likeliks_read_depth_only = matrix_of_likeliks
  
  sizesOfPointsFromLocalSds <- 0.1 / localSds 
  
  
  
  numberOfCNVsIsSufficientlySmall = F
  iterations = 0
  maxIteration = opt$maxNumIter
  while (!numberOfCNVsIsSufficientlySmall & iterations < maxIteration) {
    found_CNVs_total <- matrix(0, nrow=0, ncol=6)
    colnames(found_CNVs_total) <- c("#chr", "start", "end", "CN_change", "loglikelihood", "genes")
    
    iterations = iterations + 1
    for (l in 1:length(left_borders)) {
      chrom = names(left_borders)[l]
      start = left_borders[[l]]
      end = right_borders[[l]]
      for (k in 1:2) {
        if (nrow(found_CNVs_total) > opt$maxNumGermCNVs) {
          break
        }
        output_of_plots <-  paste0(folder_name, sample_name)
        which_to_allow <- "NA"
        if (k == 1) {
          which_to_allow = which(bedFile[,1] == chrom & bedFile[,2] <= as.numeric(start) )
        } else {
          which_to_allow = which(bedFile[,1] == chrom & bedFile[,2] >= as.numeric(end) )
        }
        toyMatrixOfLikeliks = matrix_of_likeliks[which_to_allow,]
        toyBedFile = bedFile[which_to_allow,]
        found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, 3))
        toyLogFoldChange = coverage.normalised[which_to_allow,sam_no]
        toySizesOfPointsFromLocalSds = sizesOfPointsFromLocalSds[which_to_allow]
        
        
        
        ### IGV PLOTTING
        if(opt$debug) {
          print("START OF IGV PLOTTING")
        }
        
        outputFileNameCNVs <- paste0(folder_name, sample_name, "/", sample_name, "_cnvs.seg")
        outputFileNameDots <- paste0(folder_name, sample_name, "/", sample_name, "_cov.seg")
        reverseFunctionUsedToTransform = function(x) {return((2 * x ** 2))}
        outputSegmentsAndDotsFromListOfCNVs(toyBedFile, found_CNVs, start, end, outputFileNameCNVs, 
                                            outputFileNameDots, sample_name, toyLogFoldChange, reverseFunctionUsedToTransform, cn_states)
        if(opt$debug) {
          print("END OF IGV PLOTTING")
        }
        ### END OF IGV PLOTTING
        
        
        
        if (nrow(found_CNVs) > 0) {
          # UNCOMMENT FOR PLOTTING!!!
          
          cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyLogFoldChange, toyBedFile, output_of_plots, chrom, cn_states, toySizesOfPointsFromLocalSds, plottingOfPNGs)
          if (found_CNVs[1,1] != -1000) {
            found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
            if (nrow(found_CNVs_total) > opt$maxNumGermCNVs) {
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
      if (nrow(found_CNVs_total) > opt$maxNumGermCNVs & iterations != maxIteration) {
        break
      }
    }
    if (nrow(found_CNVs_total) < opt$maxNumGermCNVs) {
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
  fileToOut <- paste0(folder_name, sample_name, "/CNAs.txt")
  fileConn<-file(fileToOut)
  writeLines(c(paste("##"," QC ", finalPValue, collapse = " ")), fileConn)
  close(fileConn)
  write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)
}

