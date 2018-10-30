#!/usr/bin/env Rscript
set.seed(100)
### PART WITH PARSING OPTIONS
library("optparse")

option_list = list(
  make_option(c("-norm", "--normal"), type="character", default=NULL, 
              help="path to table with normal coverages", metavar="character"),
  
  make_option(c("-t", "--tumor"), type="character", default=NULL, 
              help="path to table with tumor coverages", metavar="character"),
  
  make_option(c("-normOff", "--normalOfftarget"), type="character", default=NULL, 
              help="path to table with normal offtarget coverages", metavar="character"),
  
  make_option(c("-tOff", "--tumorOfftarget"), type="character", default=NULL, 
              help="path to table with tumor offtarget coverages", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="./result/", 
              help="output folder path [default= %default]", metavar="character"),
  
  make_option(c("-p", "--pair"), type="character", default="pairs.txt", 
              help="file with pairing information, 1st column = tumor, 2nd column = normal [default= %default]", metavar="character"),
  
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="bed file with panel description (chr \t start \t end \t gc_content \t annotation). has to use same notation as .cov files.", metavar="character"),
  
  make_option(c("-bOff", "--bedOfftarget"), type="character", default=NULL, 
              help="offtarget bed file with panel description (chr \t start \t end \t gc_content \t annotation). has to use same notation as .cov files.", metavar="character"),
  
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
  
  make_option(c("-ss", "--scoreS"), type="double", default="100", 
              help="minimum threshold for significance somatic variants", metavar="character"),
  
  make_option(c("-ls", "--lengthS"), type="double", default="5", 
              help="minimum threshold for length of somatic variants", metavar="character"),
  
  make_option(c("-mnaxnumg", "--maxNumGermCNVs"), type="double", default="100", 
              help="maximum number of germline CNVs allowed (increase thresholds if does not meet criteria)", metavar="character"),
  
  make_option(c("-mnaxnums", "--maxNumSomCNAs"), type="double", default="100", 
              help="maximum number of somatic CNAs allowed (increase thresholds if does not meet criteria)", metavar="character"),
  
  make_option(c("-mnaxnumit", "--maxNumIter"), type="double", default="3", 
              help="maximum number of iterations of variant calling", metavar="character"),
  
  make_option(c("-bafF", "--bafFolder"), type="character", default=NULL, 
              help="folder where you put BAF frequencies (one per normal, one per tumor sample)", metavar="character"),
  
  make_option(c("-d","--debug"), action="store_true", default=FALSE, help="Print debugging information while running.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


### TESTING PART
opt$bed = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/ssSC_v4.annotated.bed"
opt$tumor = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/tumor_ontarget_v4.cov"
opt$normal = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/ontarget_v4.cov"
opt$colNum = 4
opt$pair = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/pairsNew.txt"
opt$out = "/Users/gdemidov/Tuebingen/clinCNV_dev/results"
opt$folderWithScript = "/Users/gdemidov/Tuebingen/clinCNV_dev/ClinCNV/somatic"
opt$reanalyseCohort = F
opt$bedOfftarget = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/annotated_offtarget_v4.bed"
opt$tumorOfftarget = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/tumor_offtarget_v4.cov"
opt$normalOfftarget = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/offtarget_v4.cov"
opt$bafFolder = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/baf"

### PLOTTING OF PICTURES (DOES NOT REALLY NECESSARY IF YOU HAVE IGV SEGMENTS)
plottingOfPNGs = F

if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}

framework = "germline"
if (!is.null(opt$tumor)) {
  print("Tumor file was provided. Framework is switched to somatic.")
  framework = "somatic"
}

frameworkOff = "ontarget"
if (!is.null(opt$tumorOfftarget) & !is.null(opt$normalOfftarget) & !is.null(opt$bedOfftarget)) {
  print("Offtarget files are provided. We try to utilise off-target coverage also.")
  frameworkOff = "offtarget"
}

frameworkDataTypes = "covdepth"
if (!is.null(opt$bafFolder)) {
  print("Folder with BAFs were provided. Framework swithced to BAF.")
  frameworkDataTypes = "covdepthBAF"
}


### PART WITH LIBRARIES
library(robustbase)
library(MASS)
library("data.table")
library(foreach)
library(doParallel)

no_cores <- min(detectCores() - 1, 4)
no_cores = 4
cl<-makeCluster(no_cores)
registerDoParallel(cl)



### READING DATA
setwd(opt$folderWithScript)
bedFile <- read.table(opt$bed, stringsAsFactors = F, sep="\t")
colnames(bedFile) <- c("chr.X", "start", "end", "gc")
bedFile <- bedFile[order(bedFile$chr.X, as.numeric(bedFile$start)),]


for (i in 1:20) {
  tableOfValues <- table(round(as.numeric(as.character(bedFile[,4])) / i, digits = 2) * i)
  if(sum(tableOfValues[which(tableOfValues > 100)]) / sum(tableOfValues) > 0.95) break 
}
bedFile[,4] <- round(as.numeric(as.character(bedFile[,4])) / i, digits = 2) * i


normal <- read.table(opt$normal, header=T, stringsAsFactors = F)
normal <- normal[order(normal$X.chr, as.numeric(normal$start)),]
normal <- as.matrix(normal[,opt$colNum:ncol(normal)])

if (framework == "somatic") {
  tumor <- read.table(opt$tumor, header=T, stringsAsFactors = F)
  tumor <- tumor[order(tumor$X.chr, as.numeric(tumor$start)),]
  tumor <- as.matrix(tumor[,opt$colNum:ncol(tumor)])
}

if (frameworkOff == "offtarget") {
  bedFileOfftarget <- read.table(opt$bedOfftarget, stringsAsFactors = F, sep="\t")
  colnames(bedFileOfftarget) <- c("chr.X", "start", "end", "gc")
  bedFileOfftarget <- bedFileOfftarget[order(bedFileOfftarget$chr.X, as.numeric(bedFileOfftarget$start)),]
  
  bedFileOfftarget[,4] <- round(as.numeric(as.character(bedFileOfftarget[,4])), digits = 2)
  
  normalOff <- read.table(opt$normalOfftarget, header=T, stringsAsFactors = F)
  normalOff <- normalOff[order(normalOff$X.chr, as.numeric(normalOff$start)),]
  normalOff <- as.matrix(normalOff[,opt$colNum:ncol(normalOff)])
  # remain only samples that are in Normal cohort 
  normalOff <- normalOff[,which(colnames(normalOff) %in% colnames(normal))]
  
  tumorOff <- read.table(opt$tumorOfftarget, header=T, stringsAsFactors = F)
  tumorOff <- tumorOff[order(tumorOff$X.chr, as.numeric(tumorOff$start)),]
  tumorOff <- as.matrix(tumorOff[,opt$colNum:ncol(tumorOff)])
  # remain only samples that are in Tumor cohort 
  tumorOff <- tumorOff[,which(colnames(tumorOff) %in% colnames(tumor))]
}

setwd(opt$folderWithScript)
source("generalHelpers.R")
rowsToRemove <- cleanDatasetFromLowCoveredFiles(normal)
bedFile <- bedFile[-rowsToRemove,]
normal <- normal[-rowsToRemove,]
if (framework == "somatic")
  tumor <- tumor[-rowsToRemove,]

if (frameworkOff == "offtarget") {
  rowsToRemove <- cleanDatasetFromLowCoveredFiles(normalOff)
  normalOff = normalOff[-rowsToRemove,]
  tumorOff = tumorOff[-rowsToRemove,]
  bedFileOfftarget = bedFileOfftarget[-rowsToRemove,]
}

### GC CONTENT NORMALIZATION

pairs <- read.table(opt$pair, sep=",", stringsAsFactors = F)
pairs <- data.frame(pairs, ncol=2)
pairs <- unique(pairs)


if (frameworkDataTypes == "covdepthBAF") {
  setwd(opt$folderWithScript)
  source("bafSegmentation.R",local=TRUE)
  if (!dir.exists(file.path(opt$bafFolder, "/result"))) {
    dir.create(file.path(opt$bafFolder, "/result"))
  }
  listOfValues <- returnAllowedChromsBaf(pairs, normal, tumor, opt$bafFolder)
  allowedChromsBaf <- listOfValues[[1]]
  bAlleleFreqsAllSamples <- listOfValues[[2]]
}
setwd(opt$folderWithScript)

### ON TARGET GC NORMALIZATION
lst <- gc_and_sample_size_normalise(bedFile, normal)
normal <- lst[[1]]
if (framework == "somatic") {
  if (frameworkDataTypes == "covdepthBAF") {
    lst <- gc_and_sample_size_normalise(bedFile, tumor, allowedChroms=allowedChromsBaf)
  } else {
    lst <- gc_and_sample_size_normalise(bedFile, tumor)
  }
  tumor <- lst[[1]]
  bedFile <- lst[[2]]
} else {
  bedFile <- lst[[2]]
}

### ATTEMPT FOR PURITY AND PLODIY IMPUTATION
if (frameworkDataTypes == "covdepthBAF") {
  setwd(opt$folderWithScript)
  source("bafSegmentation.R",local=TRUE)
  purityPloidy <- returnPurityPloidy(pairs, normal, tumor, opt$bafFolder, bedFile, allowedChromsBaf)
}
setwd(opt$folderWithScript)

### OFF TARGET GC NORMALIZATION
if (frameworkOff == "offtarget") {
  lst <- gc_and_sample_size_normalise(bedFileOfftarget, normalOff)
  normalOff <- lst[[1]]
  if (frameworkDataTypes == "covdepthBAF") {
    
    lst <- gc_and_sample_size_normalise(bedFileOfftarget, tumorOff, allowedChroms=allowedChromsBaf)
  } else {
    lst <- gc_and_sample_size_normalise(bedFileOfftarget, tumorOff)
  }
  tumorOff <- lst[[1]]
  bedFileOfftarget <- lst[[2]]
}

### BED FILE OFFTARGET MAY NOT CONTAIN COLUMN WITH GENES
if (frameworkOff == "offtarget") {
  if (ncol(bedFileOfftarget)  == 4) {
    bedFileOfftarget <- cbind(bedFileOfftarget, rep(0, nrow(bedFileOfftarget)))
    colnames(bedFileOfftarget) <- colnames(bedFile)
  }
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
folder_name <- paste0(opt$out, "/normal/")
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




if (framework == "germline") quit()














bordersOfChroms <- getBordersOfChromosomes(bedFile)
### PROCESSING OF SOMATIC VARIANTS
setwd(opt$folderWithScript)
source("helpersSomatic.R",local=TRUE)

listOfValue <- formilngLogFoldChange(pairs, normal, tumor)
matrixOfLogFold <- listOfValue[[1]]
dictFromColumnToTumor <- listOfValue[[2]]


sdsOfSomaticSamples <- apply(matrixOfLogFold, 2, determineSDsOfSomaticSample)


sdsOfProbes <- sapply(1:nrow(matrixOfLogFold), function(i) {determineSDsOfSomaticProbe(matrixOfLogFold[i,], i)})

listOfVarianceAndMultiplicator <- esimtateVarianceFromSampleNoise(sdsOfSomaticSamples, 100000)
esimtatedVarianceFromSampleNoise <- listOfVarianceAndMultiplicator[[2]]
multiplicator <- listOfVarianceAndMultiplicator[[1]]

if (frameworkOff == "offtarget") {
  listOfValue <- formilngLogFoldChange(pairs, normalOff, tumorOff)
  matrixOfLogFoldOff =  listOfValue[[1]]
  dictFromColumnToTumor = listOfValue[[2]]
  sdsOfSomaticSamplesOff <- apply(matrixOfLogFoldOff, 2, determineSDsOfSomaticSample)
  
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
purity <- seq(from=10, to=101, by=5) / 100
purities <- c()
copy_numbers_used <- c()

for (pur in purity) {
  for (cn in copy_numbers) {
    # VALUES NOT CN NEUTRAL
    if (cn != 2 & (!(cn == 0 & pur <= 0.5))) {
      cn_state <- (1 - pur) * 2 + pur * cn
      cn_states <- c(cn_states, cn_state)
      purities <- c(purities, pur)
      copy_numbers_used <- c(copy_numbers_used, cn)
    }
    if (cn == 2) {
      if (pur >= 0.2) {
        cn_states <- c(cn_states, 2)
        copy_numbers_used <- c(copy_numbers_used, 2)
        purities <- c(purities, pur)
      }
    }
  }
}

cn_state[which(cn_state < 0.01)] = 0.01


final_order <- order(cn_states)
cn_states <- cn_states[final_order]
copy_numbers_used = copy_numbers_used[final_order]
purities = purities[final_order]
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

purities <- c(0, purities)
copy_numbers_used <- c(2, copy_numbers_used)
cn_states <- c(2, cn_states)

# CORRECTION OF CNS!!!
normalCoverage <- rpois(1000000, lambda=50)
tumorCoverage <- rpois(1000000, lambda=50)
sd_to_normalise = sd(log2(tumorCoverage / normalCoverage))
multipliersDueToLog <- c(1)
for (state in 2:length(cn_states)) {
  tumorCoverage <- rpois(1000000, lambda=25 * cn_states[state])
  sd_to_normalise_tumor = sd(log2(tumorCoverage / normalCoverage))
  multipliersDueToLog <- c(multipliersDueToLog, sd_to_normalise_tumor / sd_to_normalise)
}
multipliersDueToLog[which(is.nan(multipliersDueToLog))] <- max(multipliersDueToLog[which(!is.nan(multipliersDueToLog))])


startCoordOfNonInterruptedSegment = 1
shiftsOfCoverage <- c()




folder_name <- paste0(opt$out, "/somatic/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

allPotentialPurities <- unique(purities)
for (sam_no in 1:ncol(matrixOfLogFold)) {
  # To speed up reiteration, we do not want match between BAF file and bed file a lot of times
  if (frameworkDataTypes == "covdepthBAF") {
    closestBedRegions <- c()
    vectorsWithRegionCoordsFilled = F
  }
  
  sample_name <- colnames(matrixOfLogFold)[sam_no]
  print(paste("We are working on sample name:", sample_name))
  print(Sys.time())
  
  
  
  
  expectedMaxPurity = 1.0
  if (exists("purityPloidy")) {
    if (sample_name %in% names(purityPloidy)) expectedMaxPurity = min(1.0, purityPloidy[[sample_name]] + 0.05)
  }
  expectedMinPurity = max(0.1, (expectedMaxPurity - 0.1) / 2)
  expectedMinPurity = allPotentialPurities[which.min(abs(allPotentialPurities - expectedMinPurity))]
  expectedMaxPurity = allPotentialPurities[which.min(abs(allPotentialPurities - expectedMaxPurity))]
  
  
  if (!dir.exists(paste0(folder_name, sample_name)) | (opt$reanalyseCohort == T)) {
    
    dir.create(paste0(folder_name, sample_name))
    setwd(paste0(folder_name, sample_name))
    
    copyNumbersInsideExpectedPurities = F
    numberOfIterations = 0
    maxNumberOfIterations = 5    
    while(!copyNumbersInsideExpectedPurities & numberOfIterations < maxNumberOfIterations ) {
      # CLEAN FOLDER IN THE BEGINNING OF EACH ITERATION
      do.call(file.remove, list(list.files(paste0(folder_name, sample_name), full.names = TRUE)))
      
      print(paste("Current max purity:", expectedMaxPurity))
      print(paste("Current min purity:", expectedMinPurity))
      
      numberOfIterations = numberOfIterations + 1
      print(paste("Current iteration:", numberOfIterations))
      ### WE MAKE A ROUGH GUESS ON PURITY, BUT IF FOUND VARIANTS DO NOT FIT - WE INCREASE THE BORDERS
      
      maxPurityDifferent = F
      minPurityDifferent = F
      
      # BASED ON PURITY WE CAN EXCLUDE SOME COPY NUMBER
      topBorder <- 2 * (1 - expectedMinPurity) + expectedMinPurity * 3
      bottomBorder <-  2 * (1 - expectedMinPurity) + expectedMinPurity * 1
      indices = which(cn_states >= topBorder | cn_states <= bottomBorder | cn_states == 2)
      local_purities <- purities[indices]
      local_copy_numbers_used <- copy_numbers_used[indices]
      local_cn_states <- cn_states[indices]
      local_multipliersDueToLog <- multipliersDueToLog[indices]
      
      indices_to_remove_by_purity <- which((purities > expectedMaxPurity | purities < expectedMinPurity) & purities > 0)
      local_purities <- purities[-indices_to_remove_by_purity]
      local_copy_numbers_used <- copy_numbers_used[-indices_to_remove_by_purity]
      local_cn_states <- cn_states[-indices_to_remove_by_purity]
      local_multipliersDueToLog <- multipliersDueToLog[-indices_to_remove_by_purity]
      
      # PART FOR MATRIX OF CLONALITY (ONLY 2 CLONES)
      uniqueLocalPurities = unique(local_purities)
      zeroPurity <- which(uniqueLocalPurities == 0)
      uniqueLocalPurities = sort(uniqueLocalPurities[-zeroPurity])
      likeliksFoundCNVsVsPuritiesGlobal = matrix(nrow=0, ncol=length(uniqueLocalPurities))
      
      
      pvalsForQC <- c()
      threshold = opt$scoreS
      minimum_length_of_CNV = opt$lengthS
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
      
      
      if(opt$debug) {
        print(sam_no)
        print(sample_name)
      }
      
      
      
      dict_to_output = c()
      
      
      
      
      
      #### CORRECTION - IF THE SAMPLE HAS TOO MANY CNAS, WE EXPECT SOME SHIFT THERE
      if (frameworkDataTypes == "covdepthBAF") {
        sampleName2 <- strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
        position <- which(startsWith(names(allowedChromsBaf), prefix=sampleName2))
        if (length(position) == 1) {
          allowedChromsBafSample <- allowedChromsBaf[[position]]
          
          if (!sampleInOfftarget) {
            seq_to_exclude <- which(!bedFile[,1] %in% allowedChromsBafSample)
            shiftOfCoverage <- median(matrixOfLogFold[-seq_to_exclude,sam_no])
          } else {
            sam_no_off = which(colnames(matrixOfLogFoldOff) == sample_name)
            globalBed <- rbind(bedFile, bedFileOfftarget)
            globalLogFold <- c( matrixOfLogFold[,sam_no], matrixOfLogFoldOff[,sam_no_off])
            seq_to_exclude <- which(!globalBed[,1] %in% allowedChromsBafSample)
            shiftOfCoverage <- median(globalLogFold[-seq_to_exclude])
          }
          shiftsOfCoverage <- c(shiftsOfCoverage, shiftOfCoverage)
          matrixOfLogFold[,sam_no] = matrixOfLogFold[,sam_no] - shiftOfCoverage
          if (sampleInOfftarget)
            matrixOfLogFoldOff[,sam_no_off] = matrixOfLogFoldOff[,sam_no_off] - shiftOfCoverage
        }
        
      }
      
      
      
      
      
      matrix_of_likeliks <- form_matrix_of_likeliks_one_sample(1, ncol(matrixOfLogFold), matrixOfLogFold[,sam_no], localSds, log2(local_cn_states/2), local_multipliersDueToLog)
      
      
      ### ADD LIKELIHOODS
      if (frameworkDataTypes == "covdepthBAF") {
        print("Started BAF calculation")
        print(Sys.time())
        
        numberOfAssignedPositions = 0
        sampleName2 <- strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1]
        position <- which(startsWith(names(allowedChromsBaf), prefix=sampleName2))
        if (length(position) == 1) {
          bAlleleFreqsTumor <- bAlleleFreqsAllSamples[[position]][[ strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][1] ]]
          bAlleleFreqsNormal <- bAlleleFreqsAllSamples[[position]][[ strsplit(colnames(matrixOfLogFold)[sam_no], split="-")[[1]][2] ]]
          if (length(closestBedRegions) == 0) closestBedRegions = rep(0, nrow(bAlleleFreqsTumor))
          for (i in 1:nrow(bAlleleFreqsTumor)) {
            if (vectorsWithRegionCoordsFilled) {
              closestBedRegion = closestBedRegions[i]
            } else {
              closestBedRegion <- which(bedFile[,1] == bAlleleFreqsTumor[i,1] & bedFile[,2] - 250 <= bAlleleFreqsTumor[i,2] & bedFile[,3] + 250 >= bAlleleFreqsTumor[i,3])
              if (length(closestBedRegion) >= 1) {
                closestBedRegion = closestBedRegion[1]
                closestBedRegions[i] = closestBedRegion
              } else {
                # if there is no matching position we put 0
                closestBedRegions[i] = 0
                closestBedRegion = 0
              }
            }
            altAlleleDepth <- as.numeric(bAlleleFreqsTumor[i,5])
            overallDepth <- round(as.numeric(bAlleleFreqsTumor[i,6]))
            altAlleleDepth = round(altAlleleDepth * overallDepth)
            
            if (length(closestBedRegion) == 1 & closestBedRegion != 0) {
              numberOfAssignedPositions = numberOfAssignedPositions + 1
              
              vecOfLikeliks <- rep(0, ncol(matrix_of_likeliks))
              for (j in 1:length(local_cn_states)) {
                pur = local_purities[j]
                cn = local_copy_numbers_used[j]
                likelihood = -2 * likelihoodOfSNVBasedOnCN(altAlleleDepth, overallDepth, pur, cn)
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
      
      
      sizesOfPointsFromLocalSds <- 0.5 / localSds 
      if (sampleInOfftarget) {
        matrix_of_likeliks_off <- form_matrix_of_likeliks_one_sample(1, ncol(matrixOfLogFoldOff), matrixOfLogFoldOff[,sam_no_off], localSdsOff, log2(local_cn_states/2), local_multipliersDueToLog)
        globalMatrOfLikeliks <- rbind(matrix_of_likeliks, matrix_of_likeliks_off)
        globalBed <- rbind(bedFile, bedFileOfftarget)
        sizesOfPointsFromLocalSdsOff <- 0.5 / localSdsOff
        vecOfOrder = order(globalBed[,1], globalBed[,2])
        globalSizesOfPoints <- c(sizesOfPointsFromLocalSds, sizesOfPointsFromLocalSdsOff)[vecOfOrder]
        globalMatrOfLikeliks <- globalMatrOfLikeliks[vecOfOrder,]
        globalLogFold <- c( matrixOfLogFold[,sam_no], matrixOfLogFoldOff[,sam_no_off])[vecOfOrder]
        globalSds <-  c(localSds, localSdsOff)[vecOfOrder]
        globalBed <- globalBed[vecOfOrder,]
      }
      
      
      
      
      
      found_CNVs_total <- matrix(0, nrow=0, ncol=8)
      colnames(found_CNVs_total) <- c("#chr", "start", "end", "tumor_CN_change", "tumor_purity", "absolute_CN_change", "loglikelihood", "genes")
      allDetectedPurities = c()
      for (l in 1:length(left_borders)) {
        
        
        chrom = names(left_borders)[l]
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
                  likeliksFoundCNVsVsPurities[q, m] = min(apply(toyMatrixOfLikeliks[(startOfCNV + 1):(endOfCNV - 1),which(local_purities == localPurityCurrent)], 2, sum))
                }
              }
            }
            likeliksFoundCNVsVsPuritiesGlobal = rbind(likeliksFoundCNVsVsPuritiesGlobal, likeliksFoundCNVsVsPurities)
          }
          
          
          ### CHECKING OUT PURITIES OF FOUND CNVS
          if (!chrom %in% c("chrX", "chrY", "X", "Y")) {
            detectedPurities = unique(local_purities[found_CNVs[,4]])
            allDetectedPurities = unique(c(allDetectedPurities, detectedPurities))
            
            if (max(detectedPurities) == expectedMaxPurity & expectedMaxPurity != 1.0) {
              expectedMaxPurity = min(1.0, expectedMaxPurity + 0.05)
              maxPurityDifferent = T
            }
            if (min(detectedPurities) == expectedMinPurity & expectedMinPurity != 0.1) {
              expectedMinPurity = max(0.1, expectedMinPurity - 0.05)
              minPurityDifferent = T
            }
          }
          
          
          ### IGV PLOTTING
          if(opt$debug) {
            print("START OF IGV PLOTTING")
          }
          
          outputFileNameCNVs <- paste0(folder_name, sample_name, "/", sample_name, "_cnvs.seg")
          outputFileNameDots <- paste0(folder_name, sample_name, "/", sample_name, "_cov.seg")
          reverseFunctionUsedToTransform = function(x) {return((2 ** (x + 1)))}
          outputSegmentsAndDotsFromListOfCNVs(toyBedFile, found_CNVs, start, end, outputFileNameCNVs, 
                                              outputFileNameDots, sample_name, toyLogFoldChange, reverseFunctionUsedToTransform, local_cn_states)
          if(opt$debug) {
            print("END OF IGV PLOTTING")
          }
          ### END OF IGV PLOTTING
          
          
          
          if (nrow(found_CNVs) == 0 & length(which_to_allow) > 1) {
            found_CNVs = matrix(c(-1000, 1, length(which_to_allow), 1), nrow=1)
            #output_of_plots = paste0(output_of_plots, "/normal")
            #if (!dir.exists(output_of_plots)) {
            #  dir.create(output_of_plots)
            #}
          }
          
          
          
          if (nrow(found_CNVs) > 0) {
            cnvsToWriteOut <- plotFoundCNVs(found_CNVs, toyLogFoldChange, toyBedFile, output_of_plots, chrom, local_cn_states, local_copy_numbers_used, local_purities, toySizesOfPointsFromLocalSds, plottingOfPNGs)
            if (found_CNVs[1,1] != -1000) {
              found_CNVs_total = rbind(found_CNVs_total, cnvsToWriteOut)
            }
          }
        }
        
        # if after finishing some chromsome we know that BOTH purities are inaccurate - we stop
        if ((maxPurityDifferent & minPurityDifferent) & !(numberOfIterations == maxNumberOfIterations)) {
          #dir.create(paste0(folder_name, sample_name))
          #setwd(paste0(folder_name, sample_name))
          break
        }
      }
      # if after finishing all chromsomes we know that at least one purity was inaccurate - we stop
      print("All detected purities")
      print(allDetectedPurities)
      if (length(allDetectedPurities) > 0) {
        ### If we made a big mistake and OVERESTIMATED purity (no CNVs detected with such purity) - we have to reduce max purity setting to max detected one
        maxDetectedPurity = max(allDetectedPurities)
        minDetectedPurity = min(allDetectedPurities)
        print(paste("Max detected purity", maxDetectedPurity))
        if (expectedMaxPurity > maxDetectedPurity + 0.05) {
          maxPurityDifferent = T
          expectedMaxPurity = maxDetectedPurity + 0.05
          expectedMinPurity = min(max(0.1, minDetectedPurity - 0.05), max(0.1, (expectedMaxPurity - 0.1) / 2))
          expectedMinPurity = allPotentialPurities[which.min(abs(allPotentialPurities - expectedMinPurity))]
        }
        if ((maxPurityDifferent | minPurityDifferent) & !(numberOfIterations == maxNumberOfIterations)) {
          dir.create(paste0(folder_name, sample_name))
          setwd(paste0(folder_name, sample_name))
          next
        }
      }
      if ((!maxPurityDifferent & !minPurityDifferent) | numberOfIterations == maxNumberOfIterations) {
        copyNumbersInsideExpectedPurities = T
        print(paste("Converged after", numberOfIterations, "iterations. Max purity:", max(allDetectedPurities), "min purity: ", min(allDetectedPurities)))
        print("======================================")
        
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
        }
        minPointToNormalise = min(matrixOfClonality)
        matrixOfClonality[which(matrixOfClonality > minPointToNormalise + 10000)] = minPointToNormalise + 10000
        matrixOfClonality[upper.tri(matrixOfClonality)] <- NA
        hmcols<-colorRampPalette(c("blue","white","red"))(256)
        png(filename = paste0(sample_name, "_clonality.png"),
            width = 640, height = 640)
        heatmap((matrixOfClonality), scale="none", Rowv = NA, Colv = NA, col=hmcols, main=sample_name)
        dev.off()
      }
      if (length(pvalsForQC > 1)) {
        finalPValue <- 0
      } else {
        finalPValue = 0
      }
      fileToOut <- paste0(folder_name, sample_name, "/CNAs.txt")
      fileConn<-file(fileToOut)
      writeLines(c(paste("##"," QC ", 0, "purity by BAF (if != 1):", round(expectedMaxPurity - 0.05, digits=3), collapse = " ")), fileConn)
      close(fileConn)
      if(opt$debug) {
        print(found_CNVs_total)
      }
      write.table(found_CNVs_total, file = fileToOut, quote=F, row.names = F, sep="\t", append = T)	
    }  
  }
}
