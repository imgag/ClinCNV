#!/usr/bin/env Rscript
set.seed(100)
### PART WITH PARSING OPTIONS
library("optparse")

option_list = list(
  make_option(c("-norm", "--normal"), type="character", default=NULL, 
              help="path to table with normal coverages", metavar="character"),
  
  make_option(c("-t", "--tumor"), type="character", default=NULL, 
              help="path to table with tumor coverages", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="./result/", 
              help="output folder path [default= %default]", metavar="character"),
  
  make_option(c("-p", "--pair"), type="character", default="pairs.txt", 
              help="file with pairing information, 1st column = tumor, 2nd column = normal [default= %default]", metavar="character"),
  
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="bed file with panel description (chr \t start \t end \t gc_content \t annotation). has to use same notation as .cov files.", metavar="character"),
  
  make_option(c("-num", "--colNum"), type="double", default=4, 
              help="column where coverages start", metavar="character"),
  
  make_option(c("-script", "--folderWithScript"), type="character", default="./", 
              help="folder where you put script", metavar="character"),
  
  make_option(c("-r", "--reanalyseCohort"), type="logical", default=F, 
              help="if T, reanalyses whole cohort [default= %default]", metavar="character"),
  
  make_option(c("-sg", "--scoreG"), type="double", default="60", 
              help="minimum threshold for significance germline variants", metavar="character"),
  
  make_option(c("-lg", "--lengthG"), type="double", default="3", 
              help="minimum threshold for length of germline variants", metavar="character"),
  
  make_option(c("-ss", "--scoreS"), type="double", default="60", 
              help="minimum threshold for significance somatic variants", metavar="character"),
  
  make_option(c("-ls", "--lengthS"), type="double", default="5", 
              help="minimum threshold for length of somatic variants", metavar="character"),
  
  make_option(c("-mnaxnumg", "--maxNumGermCNVs"), type="double", default="100", 
              help="maximum number of germline CNVs allowed (increase thresholds if does not meet criteria)", metavar="character"),
  
  make_option(c("-mnaxnums", "--maxNumSomCNAs"), type="double", default="100", 
              help="maximum number of somatic CNAs allowed (increase thresholds if does not meet criteria)", metavar="character"),
  
  make_option(c("-mnaxnumit", "--maxNumIter"), type="double", default="3", 
              help="maximum number of iterations of variant calling", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

### TESTING PART
opt$bed = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/ssSC_v2.annotated.bed"
opt$tumor = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/tumor2.cov"
opt$normal = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/normal2.cov"
opt$colNum = 4
opt$pair = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/pairs.txt"
opt$out = "/Users/gdemidov/Tuebingen/forFranz/results/"
opt$folderWithScript = "/Users/gdemidov/Tuebingen/forFranz/ClinCNV/somatic"
opt$reanalyseCohort = T
opt$scoreS = 40
opt$lengthS = 1

if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}

framework = "germline"
if (!is.null(opt$tumor)) {
  print("Tumor file was provided. Framework is switched to somatic.")
  framework = "somatic"
}


### PART WITH LIBRARIES
library(robust)
library(robustbase)
library(MASS)
library("data.table")
library(foreach)
library(doParallel)

no_cores <- min(detectCores() - 1, 4)
no_cores = 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)



### READING DATA

setwd(opt$folderWithScript)
bedFile <- read.table(opt$bed, stringsAsFactors = F, sep="\t")
colnames(bedFile) <- c("chr.X", "start", "end", "gc")
bedFile <- bedFile[order(bedFile$chr.X, as.numeric(bedFile$start)),]

bedFile[,4] <- round(as.numeric(as.character(bedFile[,4])), digits = 2)


normal <- read.table(opt$normal, header=T, stringsAsFactors = F)
normal <- normal[order(normal$X.chr, as.numeric(normal$start)),]
normal <- as.matrix(normal[,opt$colNum:ncol(normal)])

if (framework == "somatic") {
  tumor <- read.table(opt$tumor, header=T, stringsAsFactors = F)
  tumor <- tumor[order(tumor$X.chr, as.numeric(tumor$start)),]
  tumor <- as.matrix(tumor[,opt$colNum:ncol(tumor)])
}







medians <- apply(sqrt(normal), 1, median)
rowsToRemove <- which(medians < 0.2)
bedFile <- bedFile[-rowsToRemove,]
normal <- normal[-rowsToRemove,]
if (framework == "somatic")
  tumor <- tumor[-rowsToRemove,]

### GC CONTENT NORMALIZATION
setwd(opt$folderWithScript)
source("generalHelpers.R")


lst <- gc_and_sample_size_normalise(bedFile, normal)
normal <- lst[[1]]
if (framework == "somatic") {
  lst <- gc_and_sample_size_normalise(bedFile, tumor)
  tumor <- lst[[1]]
  bedFile <- lst[[2]]
} else {
  bedFile <- lst[[2]]
}


### EXTRACTING INFORMATION FROM BED
bordersOfChroms <- getBordersOfChromosomes(bedFile)


### PROCESSING OF GERMLINE VARIANTS
setwd(opt$folderWithScript)
source("helpersGermline.R")
coverage <- sqrt(as.matrix(normal))



medians <- sapply(1:nrow(coverage), function(i) {EstimateModeSimple(coverage[i,], bedFile[i,1])})
coverage.normalised = sweep(coverage, 1, medians, FUN="/")
coverage.normalised <- coverage.normalised[, order((colnames(coverage.normalised)))]




sdsOfProbes <- sapply(1:nrow(coverage.normalised), function(i) {determineSDsOfGermlineProbe(coverage.normalised[i,], i)})

# In exome seq it is often the case that some hypervariable regions cause false positive calls.
# We remove all probes that look suspicious to us
# Moreover - probes with huge variability does not allow detection of CNVs and are useless
threshold <- min(quantile(sdsOfProbes, 0.99), 0.5)
probesToRemove <- which(sdsOfProbes > threshold)
coverage.normalised <- coverage.normalised[-probesToRemove,]
bedFile <- bedFile[-probesToRemove,]
sdsOfProbes <- sdsOfProbes[-probesToRemove]
normal <- normal[-probesToRemove,]
tumor <- tumor[-probesToRemove,]






autosomes <- which(!bedFile[,1] %in% c("chrX", "chrY", "X", "Y"))
sdsOfGermlineSamples <- apply(coverage.normalised[autosomes,], 2, determineSDsOfGermlineSample)





#locationsShiftedLogFoldChanges <- sweep(matrixOfLogFold, 1, locations)

listOfVarianceAndMultiplicator <- esimtateVarianceFromSampleNoise(sdsOfGermlineSamples, 30000)
esimtatedVarianceFromSampleNoise <- listOfVarianceAndMultiplicator[[2]]
multiplicator <- listOfVarianceAndMultiplicator[[1]]

vect_of_t_likeliks <- fast_dt_list(ncol(coverage.normalised) - 1)

cn_states <- 0:20



cytobands <- read.table("cytobands.txt",stringsAsFactors = F, header = F, sep="\t", row.names=NULL)
left_borders <- vector(mode="list", length=nrow(cytobands)/2)
right_borders <- vector(mode="list", length=nrow(cytobands)/2)
odd_numbers <- seq(from=1, to=nrow(cytobands), by=2)
even_numbers <- seq(from=2, to=nrow(cytobands), by=2)
names(left_borders) = as.vector(t(cytobands[,1]))[odd_numbers]
names(right_borders) = as.vector(t(cytobands[,1]))[even_numbers]
for (i in 1:length(odd_numbers)) {
  left_borders[[i]] = cytobands[odd_numbers[i], 2]
  right_borders[[i]] = cytobands[even_numbers[i], 3]
}




startCoordOfNonInterruptedSegment = 1
shiftsOfCoverage <- c()
colours = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414,30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414)]

overallResult <- matrix(0, nrow=0, ncol=7)
folder_name <- paste0(opt$out, "normal/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}


for (sam_no in 1:ncol(coverage.normalised)) {
  threshold = opt$scoreG
  minimum_length_of_CNV = opt$lengthG
  price_per_tile = 1
  initial_state <- 3
  
  
  localSds = sdsOfProbes * esimtatedVarianceFromSampleNoise[sam_no] * multiplicator
  
  
  
  sample_name <- colnames(coverage.normalised)[sam_no]
  print(sam_no)
  print(sample_name)
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
  
        if (nrow(found_CNVs) > 0) {
          # UNCOMMENT FOR PLOTTING!!!
          
          cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyLogFoldChange, toyBedFile, output_of_plots, chrom, cn_states, toySizesOfPointsFromLocalSds)
          if (found_CNVs[1,1] != -1000) {
            found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
              if (nrow(found_CNVs_total) > opt$maxNumGermCNVs) {
                break
              }
            }
          for (i in 1:nrow(found_CNVs)) {
  
            CNVnamesInside <- unlist(unique(toyBedFile[found_CNVs[i,2]:found_CNVs[i,3],4]))
            print(CNVnamesInside)
            
            
            CNVentry = matrix(c(sample_name, chrom, toyBedFile[found_CNVs[i,2],2], toyBedFile[found_CNVs[i,3],3], 
                                paste(CNVnamesInside, collapse=", "),
                                found_CNVs[i,4] - 1, 
                                found_CNVs[i,5]),
                              nrow=1)
            print(CNVentry)
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

































if (framework == "germline") quit()

### PROCESSING OF SOMATIC VARIANTS
setwd(opt$folderWithScript)
source("helpersSomatic.R")
pairs <- read.table(opt$pair, sep=",", stringsAsFactors = F)
pairs <- data.frame(pairs, ncol=2)
pairs <- unique(pairs)

matrixOfLogFold <- formilngLogFoldChange(pairs)


sdsOfSomaticSamples <- apply(matrixOfLogFold, 2, determineSDsOfSomaticSample)

sdsOfProbes <- sapply(1:nrow(matrixOfLogFold), function(i) {determineSDsOfSomaticProbe(matrixOfLogFold[i,], i)})


listOfVarianceAndMultiplicator <- esimtateVarianceFromSampleNoise(sdsOfSomaticSamples, 100000)
esimtatedVarianceFromSampleNoise <- listOfVarianceAndMultiplicator[[2]]
multiplicator <- listOfVarianceAndMultiplicator[[1]]

### FORMING MATRIX OF LIKELIHOODS
vect_of_t_likeliks <- fast_dt_list(ncol(matrixOfLogFold) - 1)




cn_states <- c()
purity <- seq(from=10, to=101, by=5) / 100
for (pur in purity) {
  cn_state <- 0 + 1 * pur * seq(from=0, to=20, by=1)
  cn_states <- c(cn_states, cn_state)
}
cn_states <- unique(cn_states[-which(cn_states > 1.5 & cn_states < 2.8)])
cn_state[which(cn_state < 0.1)] = 0.1

cn_states <- sort(cn_states)
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

cn_states <- c(2, cn_states)


# CORRECTION OF CNS!!!
normalCoverage <- rpois(100000, lambda=50)
tumorCoverage <- rpois(100000, lambda=50)
sd_to_normalise = sd(log2(tumorCoverage / normalCoverage))
multipliersDueToLog <- c(1)
for (state in 2:length(cn_states)) {
  tumorCoverage <- rpois(100000, lambda=50 * cn_states[state])
  sd_to_normalise_tumor = sd(log2(tumorCoverage / normalCoverage))
  multipliersDueToLog <- c(multipliersDueToLog, sd_to_normalise_tumor / sd_to_normalise)
}
multipliersDueToLog[which(is.nan(multipliersDueToLog))] <- max(multipliersDueToLog[which(!is.nan(multipliersDueToLog))])







startCoordOfNonInterruptedSegment = 1
shiftsOfCoverage <- c()




folder_name <- paste0(opt$out, "somatic/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

deletedSamples <- c()
for (sam_no in 1:ncol(matrixOfLogFold)) {
  pvalsForQC <- c()
  threshold = opt$scoreS
  minimum_length_of_CNV = opt$lengthS
  price_per_tile = 0.1
  initial_state <- 1
  
  
  localSds = sdsOfProbes * esimtatedVarianceFromSampleNoise[sam_no] * multiplicator
  
  
  
  sample_name <- colnames(matrixOfLogFold)[sam_no]
  print(sam_no)
  print(sample_name)
  if (!dir.exists(paste0(folder_name, sample_name)) | (opt$reanalyseCohort == T)) {
    dir.create(paste0(folder_name, sample_name))
    
    setwd(paste0(folder_name, sample_name))
    
    
    dict_to_output = c()
    matrix_of_likeliks <- form_matrix_of_likeliks_one_sample(1, ncol(matrixOfLogFold), sam_no, localSds, matrixOfLogFold, log2(cn_states/2), multipliersDueToLog)
    
    #### CORRECTION - IF THE SAMPLE HAS TOO MANY CNAS, WE EXPECT SOME SHIFT THERE
    found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, matrix_of_likeliks, 1))
    seq_to_exclude <- c()
    if (nrow(found_CNVs) > 0)
      for (i in 1:nrow(found_CNVs)) {
        seq_to_exclude <- c(seq_to_exclude, (found_CNVs[i,2]):(found_CNVs[i,3]))
      }
    if (length(seq_to_exclude) > 0) {
      shiftOfCoverage <- median(matrixOfLogFold[-seq_to_exclude,sam_no])
      shiftsOfCoverage <- c(shiftsOfCoverage, shiftOfCoverage)
      matrixOfLogFold[,sam_no] = matrixOfLogFold[,sam_no] - shiftOfCoverage
    }
    if (sam_no %in% deletedSamples) {
      matrixOfLogFold[,sam_no] = matrixOfLogFold[,sam_no] - log2(2)
    }
    
    
    
    
    
    matrix_of_likeliks <- form_matrix_of_likeliks_one_sample(1, ncol(matrixOfLogFold), sam_no, localSds, matrixOfLogFold, log2(cn_states/2), multipliersDueToLog)
    matrix_of_likeliks_read_depth_only = matrix_of_likeliks
    
    sizesOfPointsFromLocalSds <- 0.5 / localSds 
    

    found_CNVs_total <- matrix(0, nrow=0, ncol=6)
    colnames(found_CNVs_total) <- c("#chr", "start", "end", "CN_change", "loglikelihood", "genes")
    for (l in 1:length(left_borders)) {
      
      chrom = names(left_borders)[l]
      start = left_borders[[l]]
      end = right_borders[[l]]
      for (k in 1:2) {
        output_of_plots <-  paste0(folder_name, sample_name)
        which_to_allow <- "NA"
        if (k == 1) {
          which_to_allow = which(bedFile[,1] == chrom & bedFile[,2] <= start )
        } else {
          which_to_allow = which(bedFile[,1] == chrom & bedFile[,2] >= end )
        }
        
        toyMatrixOfLikeliks = matrix_of_likeliks[which_to_allow,]
        toyBedFile = bedFile[which_to_allow,]
        found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, toyMatrixOfLikeliks, 1))
        print(found_CNVs)
        print(l)
        print(k)
        toyLogFoldChange = matrixOfLogFold[which_to_allow,]
        toySds <- localSds[which_to_allow]
        toyMultipliersDueToLog <- multipliersDueToLog[which_to_allow]
        pvalueForThisArmQC <- qcControl(sam_no, toyLogFoldChange, toySds, toyMultipliersDueToLog, found_CNVs, 0.8)
        
        if (pvalueForThisArmQC > 0) {
          pvalsForQC <- c(pvalsForQC, pvalueForThisArmQC)
        }
        
        toyLogFoldChange = matrixOfLogFold[which_to_allow,sam_no]
        toySizesOfPointsFromLocalSds = sizesOfPointsFromLocalSds[which_to_allow]
        if (nrow(found_CNVs) == 0 & length(which_to_allow) > 1) {
          found_CNVs = matrix(c(-1000, 1, length(which_to_allow), 1), nrow=1)
          output_of_plots = paste0(output_of_plots, "/normal")
          if (!dir.exists(output_of_plots)) {
            dir.create(output_of_plots)
          }
        }
        if (nrow(found_CNVs) > 0) {
          cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyLogFoldChange, toyBedFile, output_of_plots, chrom, cn_states, toySizesOfPointsFromLocalSds)
          if (found_CNVs[1,1] != -1000) {
            found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
          }
        }
      }
    }
  if (length(pvalsForQC > 1)) {
    finalPValue <- median(pvalsForQC)
  } else {
    finalPValue = "UNKNOWN"
  }
  fileToOut <- paste0(folder_name, sample_name, "/CNAs.txt")
  fileConn<-file(fileToOut)
  writeLines(c(paste("##"," QC ", finalPValue, collapse = " ")), fileConn)
  close(fileConn)
  write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)	
	
	
  }  
}
