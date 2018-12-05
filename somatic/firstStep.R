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
  
  make_option(c("-sg", "--scoreG"), type="double", default="40", 
              help="minimum threshold for significance germline variants", metavar="character"),
  
  make_option(c("-lg", "--lengthG"), type="double", default="2", 
              help="minimum threshold for length of germline variants", metavar="character"),
  
  make_option(c("-ss", "--scoreS"), type="double", default="200", 
              help="minimum threshold for significance somatic variants", metavar="character"),
  
  make_option(c("-ls", "--lengthS"), type="double", default="4", 
              help="minimum threshold for length of somatic variants", metavar="character"),
  
  make_option(c("-mnaxnumg", "--maxNumGermCNVs"), type="double", default="100", 
              help="maximum number of germline CNVs allowed (increase thresholds if does not meet criteria)", metavar="character"),
  
  make_option(c("-mnaxnums", "--maxNumSomCNAs"), type="double", default="100", 
              help="maximum number of somatic CNAs allowed (increase thresholds if does not meet criteria)", metavar="character"),
  
  make_option(c("-mnaxnumit", "--maxNumIter"), type="double", default="3", 
              help="maximum number of iterations of variant calling", metavar="character"),
  
  make_option(c("-bafF", "--bafFolder"), type="character", default=NULL, 
              help="folder where you put BAF frequencies (one per normal, one per tumor sample)", metavar="character"),
  
  make_option(c("-normS", "--normalSample"), type="character", default=NULL, 
              help="name of normal sample to analyse (if only one sample has to be analysed)", metavar="character"),
  
  make_option(c("-tumorS", "--tumorSample"), type="character", default=NULL, 
              help="name of tumor sample to analyse (if only one sample has to be analysed, normal has to be provided too)", metavar="character"),
  
  make_option(c("-d","--debug"), action="store_true", default=FALSE, help="Print debugging information while running.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);




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
library(mclust)

no_cores <- min(detectCores() - 1, 4)
no_cores = 4
cl<-makeCluster(no_cores)
registerDoParallel(cl)



### READING DATA
setwd(opt$folderWithScript)
bedFile <- read.table(opt$bed, stringsAsFactors = F, sep="\t", comment.char="&", header=F)
if (!startsWith(bedFile[,1], "chr"))
  bedFile[,1] <- paste0("chr", bedFile[,1])
colnames(bedFile) <- c("chr.X", "start", "end", "gc")
bedFile <- bedFile[order(bedFile$chr.X, as.numeric(bedFile$start)),]


for (i in 1:20) {
  tableOfValues <- table(round(as.numeric(as.character(bedFile[,4])) / i, digits = 2) * i)
  if(sum(tableOfValues[which(tableOfValues > 100)]) / sum(tableOfValues) > 0.95) break 
}
bedFile[,4] <- round(as.numeric(as.character(bedFile[,4])) / i, digits = 2) * i
whichBedIsNA <- which(is.na(bedFile[,4]) | bedFile[,3] - bedFile[,2] < 80)
if (length(whichBedIsNA) > 0)
  bedFile = bedFile[-whichBedIsNA,]

normal <- read.table(opt$normal, header=T, stringsAsFactors = F, comment.char="&" )
if (!startsWith(normal[,1], "chr"))
  normal[,1] <- paste0("chr", normal[,1])
normal <- normal[order(normal[,1], as.numeric(normal[,2])),]
normal <- as.matrix(normal[,opt$colNum:ncol(normal)])
if (length(whichBedIsNA) > 0)
  normal = normal[-whichBedIsNA,]

if (framework == "somatic") {
  tumor <- read.table(opt$tumor, header=T, stringsAsFactors = F, comment.char="&" )
  if (!startsWith(tumor[,1], "chr"))
    tumor[,1] <- paste0("chr", tumor[,1])
  tumor <- tumor[order(tumor$X.chr, as.numeric(tumor$start)),]
  tumor <- as.matrix(tumor[,opt$colNum:ncol(tumor)])
  if (length(whichBedIsNA) > 0)
    tumor = tumor[-whichBedIsNA,]
}

if (frameworkOff == "offtarget") {
  bedFileOfftarget <- read.table(opt$bedOfftarget, stringsAsFactors = F, sep="\t")
  if (!startsWith(bedFileOfftarget[,1], "chr"))
    bedFileOfftarget[,1] <- paste0("chr", bedFileOfftarget[,1])
  colnames(bedFileOfftarget) <- c("chr.X", "start", "end", "gc")
  bedFileOfftarget <- bedFileOfftarget[order(bedFileOfftarget$chr.X, as.numeric(bedFileOfftarget$start)),]
  
  for (i in 1:20) {
    tableOfValues <- table(round(as.numeric(as.character(bedFileOfftarget[,4])) / i, digits = 2) * i)
    if(sum(tableOfValues[which(tableOfValues > 100)]) / sum(tableOfValues) > 0.95) break 
  }
  bedFileOfftarget[,4] <- round(as.numeric(as.character(bedFileOfftarget[,4])) / i, digits = 2) * i
  
  
  normalOff <- read.table(opt$normalOfftarget, header=T, stringsAsFactors = F, comment.char="&" )
  if (!startsWith(normalOff[,1], "chr"))
    normalOff[,1] <- paste0("chr", normalOff[,1])
  normalOff <- normalOff[order(normalOff$X.chr, as.numeric(normalOff$start)),]
  normalOff <- as.matrix(normalOff[,opt$colNum:ncol(normalOff)])
  # remain only samples that are in Normal cohort 
  normalOff <- normalOff[,which(colnames(normalOff) %in% colnames(normal))]
  
  tumorOff <- read.table(opt$tumorOfftarget, header=T, stringsAsFactors = F, comment.char="&" )
  if (!startsWith(tumorOff[,1], "chr"))
    tumorOff[,1] <- paste0("chr", tumorOff[,1])
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

if (framework == "somatic") {
  pairs <- read.table(opt$pair, sep=",", stringsAsFactors = F)
  pairs <- data.frame(pairs, ncol=2)
  pairs <- unique(pairs)
}


## CHECK INPUT VALIDITY
if (!is.null(opt$normalSample)) {
  stopifnot(opt$normalSample %in% colnames(normal))
}
if (!is.null(opt$tumorSample)) {
  stopifnot(!is.null(opt$normalSample))
  stopifnot(opt$tumorSample %in% colnames(tumor))
  stopifnot(opt$tumorSample %in% pairs[,1])
  stopifnot(opt$normalSample %in% pairs[,2])
  coordOfNormalInPairs = which(pairs[,2] == opt$normalSample)
  stopifnot(opt$tumorSample %in% pairs[coordOfNormalInPairs,1])
}

lstOfChromBorders <- getCytobands("cytobands.txt")
left_borders <- lstOfChromBorders[[1]]
right_borders <- lstOfChromBorders[[2]]
ends_of_chroms <- lstOfChromBorders[[3]]

if (frameworkDataTypes == "covdepthBAF") {
  setwd(opt$folderWithScript)
  source("bafSegmentation.R",local=TRUE)
  if (!dir.exists(file.path(opt$bafFolder, "/result"))) {
    dir.create(file.path(opt$bafFolder, "/result"))
  }
  if (!is.null(opt$normalSample) & !is.null(opt$tumorSample)) {
    coordOfNormalInPairs = which(pairs[,2] == opt$normalSample & pairs[,1] == opt$tumorSample)
    pairsForBAF = pairs[coordOfNormalInPairs,,drop=F]
  } else {
    pairsForBAF = pairs
  }
  listOfValues <- returnAllowedChromsBaf(pairsForBAF, normal, tumor, opt$bafFolder, bedFile, left_borders, right_borders, ends_of_chroms)
  allowedChromsBaf <- listOfValues[[1]]
  bAlleleFreqsAllSamples <- listOfValues[[2]]
}
setwd(opt$folderWithScript)

### ON TARGET GC NORMALIZATION
if (max(bedFile[,3] - bedFile[,2]) / min(bedFile[,3] - bedFile[,2]) > 16)
  normal <- lengthBasedNormalization(normal, bedFile)
lst <- gc_and_sample_size_normalise(bedFile, normal)
normal <- lst[[1]]
if (framework == "somatic") {
  if (frameworkDataTypes == "covdepthBAF") {
    if (max(bedFile[,3] - bedFile[,2]) / min(bedFile[,3] - bedFile[,2]) > 16)
      tumor <- lengthBasedNormalization(tumor, bedFile, allowedChroms=allowedChromsBaf)
    lst <- gc_and_sample_size_normalise(bedFile, tumor, allowedChroms=allowedChromsBaf)
  } else {
    tumor <- lengthBasedNormalization(tumor, bedFile)
    lst <- gc_and_sample_size_normalise(bedFile, tumor)
  }
  tumor <- lst[[1]]
  bedFile <- lst[[2]]
} else {
  bedFile <- lst[[2]]
}



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

# FILTER LOW COVERED REGIONS
regionsToFilerOutOn <- c()
for (i in 1:nrow(normal)) {
if (framework == "somatic"){
  if (min(median(normal[i,], tumor[i,])) < 0.3) {
    regionsToFilerOutOn <- c(regionsToFilerOutOn, i)
  }
 }
}
if (length(regionsToFilerOutOn)>0)
{
	normal = normal[-regionsToFilerOutOn,] + 10**-20
	tumor = tumor[-regionsToFilerOutOn,] + 10**-20
	bedFile = bedFile[-regionsToFilerOutOn,]
}
if (frameworkOff == "offtarget") {
  regionsToFilerOutOff <- c()
  for (i in 1:nrow(normalOff)) {
    if (min(median(normalOff[i,], tumorOff[i,])) < 0.3) {
      regionsToFilerOutOff <- c(regionsToFilerOutOff, i)
    }
  }
  normalOff = normalOff[-regionsToFilerOutOff,] + 10**-20
  tumorOff = tumorOff[-regionsToFilerOutOff,] + 10**-20
  bedFileOfftarget = bedFileOfftarget[-regionsToFilerOutOff,]
}

### EXTRACTING INFORMATION FROM BED
bordersOfChroms <- getBordersOfChromosomes(bedFile)


### PROCESSING OF GERMLINE VARIANTS
setwd(opt$folderWithScript)
source("helpersGermline.R")
coverage <- sqrt(as.matrix(normal))



medians <- sapply(1:nrow(coverage), function(i) {EstimateModeSimple(coverage[i,], bedFile[i,1])})
whichMediansAreSmall <- which(medians < 0.5)
if (length(whichMediansAreSmall) > 0) {
  coverage <- coverage[-whichMediansAreSmall,]
  bedFile <- bedFile[-whichMediansAreSmall,]
  medians <- medians[-whichMediansAreSmall]
}
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
if (framework=="somatic")
  tumor <- tumor[-probesToRemove,]


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




if (framework == "germline") quit()













### PROCESSING OF SOMATIC VARIANTS
setwd(opt$folderWithScript)
source("helpersSomatic.R",local=TRUE)

listOfValue <- formilngLogFoldChange(pairs, normal, tumor)
matrixOfLogFold <- listOfValue[[1]]
dictFromColumnToTumor <- listOfValue[[2]]

bordersOfChroms <- getBordersOfChromosomes(bedFile)
sdsOfSomaticSamples <- apply(matrixOfLogFold, 2, determineSDsOfSomaticSample)


sdsOfProbes <- sapply(1:nrow(matrixOfLogFold), function(i) {determineSDsOfSomaticProbe(matrixOfLogFold[i,], i)})

listOfVarianceAndMultiplicator <- esimtateVarianceFromSampleNoise(sdsOfSomaticSamples, 100000)
esimtatedVarianceFromSampleNoise <- listOfVarianceAndMultiplicator[[2]]
multiplicator <- listOfVarianceAndMultiplicator[[1]]

if (frameworkOff == "offtarget") {
  listOfValue <- formilngLogFoldChange(pairs, normalOff, tumorOff)
  matrixOfLogFoldOff =  listOfValue[[1]]
  dictFromColumnToTumor = listOfValue[[2]]
  bordersOfChroms <- getBordersOfChromosomes(bedFileOfftarget)
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
purities <- c(0, purities)
copy_numbers_used <- c(2, copy_numbers_used)
cn_states <- c(2, cn_states)
statesUsed <- c("normal", statesUsed)






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
sdsOfNormals <- apply(normal, 1, sd)
medianBaselineSD <- median(sdsOfNormals)

sdsPois <- c()
for (i in 1:1000) {
  samplePois <- rpois(lambda=i, n=10000)
  sdsPois <- c(sdsPois, sd(samplePois / i))
}
diffs <- abs(sdsPois - medianBaselineSD)
lambdaFromSimulation <- which.min(diffs)

normalCoverage <- rpois(1000000, lambda=lambdaFromSimulation)
tumorCoverage <- rpois(1000000, lambda=lambdaFromSimulation)
sd_to_normalise = sd(log2(tumorCoverage / normalCoverage))
multipliersDueToLog <- c(1)
for (state in 2:length(cn_states)) {
  tumorCoverage <- rpois(1000000, lambda=0.5 * lambdaFromSimulation * cn_states[state])
  sd_to_normalise_tumor = sd(log2(tumorCoverage / normalCoverage))
  multipliersDueToLog <- c(multipliersDueToLog, sd_to_normalise_tumor / sd_to_normalise)
}
multipliersDueToLog[which(is.nan(multipliersDueToLog))] <- max(multipliersDueToLog[which(!is.nan(multipliersDueToLog))])

### Make SDs for states equal
for (state in unique(cn_states)) {
  whichStatesAre = which(cn_states == state)
  if (state != 2) {
    multipliersDueToLog[whichStatesAre] = median(multipliersDueToLog[whichStatesAre])
  } else {
    multipliersDueToLog[whichStatesAre] = 1
  }
}

startCoordOfNonInterruptedSegment = 1
shiftsOfCoverage <- c()




folder_name <- paste0(opt$out, "/somatic/")
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

allPotentialPurities <- unique(purities)
for (sam_no in 1:ncol(matrixOfLogFold)) {
  sample_name <- colnames(matrixOfLogFold)[sam_no]
  
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
      
      if (finalIteration ) {
        if ((abs(max(matrixOfClonality) - min(matrixOfClonality)) > 1000)) {
          clonalSignificanceThreshold = 1000
          indices = which(matrixOfClonality == min(matrixOfClonality), arr.ind = TRUE)
          oneClone = min(matrixOfClonality[indices[1], indices[1]], matrixOfClonality[indices[2], indices[2]])
          if (abs(oneClone - matrixOfClonality[indices[1,1], indices[1,2]]) > clonalSignificanceThreshold) {
            clonalBestPurities <- rownames(which(matrixOfClonality == min(matrixOfClonality), arr.ind = TRUE))
            clonalBestPurities = as.numeric(clonalBestPurities)
          } else {
            clonalBestPurities <- rownames(which(matrixOfClonality == oneClone, arr.ind = TRUE))
            clonalBestPurities = as.numeric(clonalBestPurities)
            tableOfBestPurities <- table(rownames(which(matrixOfClonality == oneClone, arr.ind = TRUE)))
            clonalBestPurities = as.numeric(names(tableOfBestPurities[which.max(tableOfBestPurities)]))
          }
          
          clonalBestPurities <- c(as.numeric(clonalBestPurities), 0)
          indices_to_remove_by_purity <- which(!(purities %in% clonalBestPurities))
          local_purities <- purities[-indices_to_remove_by_purity]
          local_copy_numbers_used <- copy_numbers_used[-indices_to_remove_by_purity]
          local_cn_states <- cn_states[-indices_to_remove_by_purity]
          local_multipliersDueToLog <- multipliersDueToLog[-indices_to_remove_by_purity]
          local_cnv_states = local_cnv_states[-indices_to_remove_by_purity]
        }
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
            clusteredResult <- densityMclust(smoothedLogFold)
            print("Mclust finished")
            bigClusters <- which(clusteredResult$parameters$pro > 0.3)
            if (length(bigClusters) == 0) {
              shiftOfCoverage <- median(globalLogFold[allowedChromosomesAutosomesOnly])
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
            clusteredResult <- densityMclust(smoothedLogFold)
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
      
      
      
      
      
      matrix_of_likeliks <- form_matrix_of_likeliks_one_sample(1, ncol(matrixOfLogFold), matrixOfLogFold[,sam_no], localSds, log2(local_cn_states/2), local_multipliersDueToLog)
      
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
      
      
      
      
      
      found_CNVs_total <- matrix(0, nrow=0, ncol=9)
      colnames(found_CNVs_total) <- c("#chr", "start", "end", "tumor_CN_change", "tumor_clonality", "CN_change", "loglikelihood", "state", "genes")
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
          
          
          
          ### IGV PLOTTING
          if(opt$debug) {
            print("START OF IGV PLOTTING")
          }
          if (finalIteration) {
          outputFileNameCNVs <- paste0(folder_name, sample_name, "/", sample_name, "_cnvs.seg")
          outputFileNameDots <- paste0(folder_name, sample_name, "/", sample_name, "_cov.seg")
          reverseFunctionUsedToTransform = function(x) {return((2 ** (x + 1)))}
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


