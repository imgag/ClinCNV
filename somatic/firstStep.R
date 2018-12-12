#!/usr/bin/env Rscript
set.seed(100)

## CHECK R VERSION
if (!(as.numeric(version$major) >= 3 & as.numeric(version$minor) > 2.0)) {
 print("Your R version is too old. We can not guarantee stable work.")
 print(version)
}
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
  
  make_option(c("-script", "--folderWithScript"), type="character", default=getwd(), 
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
  
  make_option(c("-triosFile", "--triosFile"), type="character", default=NULL, 
              help="file with information about trios, child-father-mother", metavar="character"),
  
  make_option(c("-d","--debug"), action="store_true", default=FALSE, help="Print debugging information while running.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


setwd(opt$folderWithScript)
source("generalHelpers.R")

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
cl<-makeCluster(no_cores, type="FORK")
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
normal = normal
colnames(normal) = cutX(colnames(normal))
if (!startsWith(normal[,1], "chr"))
  normal[,1] <- paste0("chr", normal[,1])
normal <- normal[order(normal[,1], as.numeric(normal[,2])),]
normal <- as.matrix(normal[,opt$colNum:ncol(normal)])
if (length(whichBedIsNA) > 0)
  normal = normal[-whichBedIsNA,]

if (framework == "somatic") {
  tumor <- read.table(opt$tumor, header=T, stringsAsFactors = F, comment.char="&" )
  tumor = tumor
  colnames(tumor) = cutX(colnames(tumor))
  if (!startsWith(tumor[,1], "chr"))
    tumor[,1] <- paste0("chr", tumor[,1])
  tumor <- tumor[order(tumor[,1], as.numeric(tumor[,2])),]
  tumor <- as.matrix(tumor[,opt$colNum:ncol(tumor)])
  if (length(whichBedIsNA) > 0)
    tumor = tumor[-whichBedIsNA,]
}

if (frameworkOff == "offtarget") {
  bedFileOfftarget <- read.table(opt$bedOfftarget, stringsAsFactors = F, sep="\t")
  if (!startsWith(bedFileOfftarget[,1], "chr"))
    bedFileOfftarget[,1] <- paste0("chr", bedFileOfftarget[,1])
  colnames(bedFileOfftarget) <- c("chr.X", "start", "end", "gc")
  bedFileOfftarget <- bedFileOfftarget[order(bedFileOfftarget[,1], as.numeric(bedFileOfftarget[,2])),]
  
  for (i in 1:20) {
    tableOfValues <- table(round(as.numeric(as.character(bedFileOfftarget[,4])) / i, digits = 2) * i)
    if(sum(tableOfValues[which(tableOfValues > 100)]) / sum(tableOfValues) > 0.95) break 
  }
  bedFileOfftarget[,4] <- round(as.numeric(as.character(bedFileOfftarget[,4])) / i, digits = 2) * i
  
  
  normalOff <- read.table(opt$normalOfftarget, header=T, stringsAsFactors = F, comment.char="&" )
  colnames(normalOff) = cutX(colnames(normalOff))
  if (!startsWith(normalOff[,1], "chr"))
    normalOff[,1] <- paste0("chr", normalOff[,1])
  normalOff <- normalOff[order(normalOff[,1], as.numeric(normalOff[,2])),]
  normalOff <- as.matrix(normalOff[,opt$colNum:ncol(normalOff)])
  # remain only samples that are in Normal cohort 
  normalOff <- normalOff[,which(colnames(normalOff) %in% colnames(normal))]
  
  tumorOff <- read.table(opt$tumorOfftarget, header=T, stringsAsFactors = F, comment.char="&" )
  colnames(tumorOff) = cutX(colnames(tumorOff))
  if (!startsWith(tumorOff[,1], "chr"))
    tumorOff[,1] <- paste0("chr", tumorOff[,1])
  tumorOff <- tumorOff[order(tumorOff[,1], as.numeric(tumorOff[,2])),]
  tumorOff <- as.matrix(tumorOff[,opt$colNum:ncol(tumorOff)])
  # remain only samples that are in Tumor cohort 
  tumorOff <- tumorOff[,which(colnames(tumorOff) %in% colnames(tumor))]
}


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

# checkSignalToNoise <- function(matr) {
#   mediansWithoutNorm <- apply(matr, 2, median)
#   sdsWithoutNorm <- apply(matr, 2, Sn)
#   return(mediansWithoutNorm / sdsWithoutNorm)
# }
# snNormWith = checkSignalToNoise(normal)
# snTumorWith = checkSignalToNoise(tumor)
# snNormWithout = checkSignalToNoise(normal)
# snTumorWithout = checkSignalToNoise(tumor)

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
if (ncol(bedFile)  == 4) {
  bedFile <- cbind(bedFile, rep(0, nrow(bedFile)))
  colnames(bedFile) <- colnames(bedFile)
}
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

print("Processing of germline variants started.")
coverage <- sqrt(as.matrix(normal))

print("Gender determination started")
genderOfSamples <- Determine.gender(coverage, bedFile)
print("Gender succesfully determined. Plot is written in your results directory.")


clusterExport(cl, c('EstimateModeSimple', 'bedFile', 'genderOfSamples', 'coverage', "lehmanHodges", 'Qn'))
medians <- parSapply(cl=cl, 1:nrow(coverage), function(i) {EstimateModeSimple(coverage[i,], bedFile[i,1], genderOfSamples)})
whichMediansAreSmall <- which(medians < 0.5)
if (length(whichMediansAreSmall) > 0) {
  coverage <- coverage[-whichMediansAreSmall,]
  bedFile <- bedFile[-whichMediansAreSmall,]
  medians <- medians[-whichMediansAreSmall]
}
coverage.normalised = sweep(coverage, 1, medians, FUN="/")
coverage.normalised <- coverage.normalised[, order((colnames(coverage.normalised)))]


clusterExport(cl, c('coverage.normalised', 'determineSDsOfGermlineProbe'))
sdsOfProbes <- parSapply(cl=cl, 1:nrow(coverage.normalised), function(i) {determineSDsOfGermlineProbe(coverage.normalised[i,], i)})

# In exome seq it is often the case that some hypervariable regions cause false positive calls.
# We remove all probes that look suspicious to us
# Moreover - probes with huge variability does not allow detection of CNVs and are useless
threshold <- min(quantile(sdsOfProbes, 0.99), 0.5)
probesToRemove <- which(sdsOfProbes > threshold)
coverage.normalised <- coverage.normalised[-probesToRemove,]
bedFile <- bedFile[-probesToRemove,]
sdsOfProbes <- sdsOfProbes[-probesToRemove]
normal <- normal[-probesToRemove,]





autosomes <- which(!bedFile[,1] %in% c("chrX", "chrY", "X", "Y"))
sdsOfGermlineSamples <- apply(coverage.normalised[autosomes,], 2, determineSDsOfGermlineSample)





#locationsShiftedLogFoldChanges <- sweep(matrixOfLogFold, 1, locations)

listOfVarianceAndMultiplicator <- esimtateVarianceFromSampleNoise(sdsOfGermlineSamples, 100000)
esimtatedVarianceFromSampleNoise <- listOfVarianceAndMultiplicator[[2]]
multiplicator <- listOfVarianceAndMultiplicator[[1]]

vect_of_t_likeliks <- fast_dt_list(ncol(coverage.normalised) - 1)
vect_of_norm_likeliks <- fast_dnorm_list()





setwd(opt$folderWithScript)
if (!is.null(opt$triosFile)) {
  source("germlineTrioSolver.R",local=TRUE)
} else {
  source("germlineSolver.R",local=TRUE)
}
stopCluster(cl)
if (framework == "germline" | !is.null(opt$triosFile)) quit()









no_cores <- min(detectCores() - 1, 4)
no_cores = 4
cl<-makeCluster(no_cores, type="FORK")
registerDoParallel(cl)


if (framework=="somatic")
  tumor <- tumor[-probesToRemove,]



setwd(opt$folderWithScript)
source("somaticSolver.R",local=TRUE)
stopCluster(cl)