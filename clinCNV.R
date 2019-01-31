#!/usr/bin/env Rscript
set.seed(100)
options(warn=-1)

## CHECK R VERSION
if (!(as.numeric(version$major) >= 3 & as.numeric(version$minor) > 2.0)) {
 print("Your R version is too old. We can not guarantee stable work.")
 print(version)
}

### PART WITH PARSING OPTIONS
library("optparse")
library(robustbase)
library(MASS)
library("data.table")
library(foreach)
library(doParallel)
library(msm)



initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)


## DETERMINE THE PATH TO THE SCRIPT AUTOMATICALLY
current_working_dir <- script.basename


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
  
  make_option(c("-num", "--colNum"), type="integer", default=4, 
              help="column where coverages start", metavar="number"),
  
  make_option(c("-script", "--folderWithScript"), type="character", default=current_working_dir, 
              help="folder where you put script", metavar="character"),
  
  make_option(c("-r", "--reanalyseCohort"), action="store_false", 
              help="if specified, reanalyses whole cohort [default= %default]"),
  
  make_option(c("-sg", "--scoreG"), type="double", default="40", 
              help="minimum threshold for significance germline variants", metavar="number"),
  
  make_option(c("-lg", "--lengthG"), type="integer", default="2", 
              help="minimum threshold for length of germline variants", metavar="number"),
  
  make_option(c("-ss", "--scoreS"), type="double", default="200", 
              help="minimum threshold for significance somatic variants", metavar="number"),
  
  make_option(c("-ls", "--lengthS"), type="integer", default="4", 
              help="minimum threshold for length of somatic variants", metavar="number"),
  
  make_option(c("-mnaxnumg", "--maxNumGermCNVs"), type="integer", default="100", 
              help="maximum number of germline CNVs allowed (increase thresholds if does not meet criteria)", metavar="number"),
  
  make_option(c("-mnaxnums", "--maxNumSomCNAs"), type="integer", default="100", 
              help="maximum number of somatic CNAs allowed (increase thresholds if does not meet criteria)", metavar="number"),
  
  make_option(c("-mnaxnumit", "--maxNumIter"), type="integer", default=3, 
              help="maximum number of iterations of variant calling", metavar="number"),
  
  make_option(c("-bafF", "--bafFolder"), type="character", default=NULL, 
              help="folder where you put BAF frequencies (one per normal, one per tumor sample)", metavar="character"),
  
  make_option(c("-normS", "--normalSample"), type="character", default=NULL, 
              help="name of normal sample to analyse (if only one sample has to be analysed)", metavar="character"),
  
  make_option(c("-tumorS", "--tumorSample"), type="character", default=NULL, 
              help="name of tumor sample to analyse (if only one sample has to be analysed, normal has to be provided too)", metavar="character"),
  
  make_option(c("-triosFile", "--triosFile"), type="character", default=NULL, 
              help="file with information about trios, child-father-mother", metavar="character"),
  
  make_option(c("-fdrG", "--fdrGermline"), type="integer", default=0, 
              help="number of iterations for FDR check (more - better, but slower, 0 = no FDR correction)", metavar="number"),

  make_option(c("-numT", "--numberOfThreads"), type="character", default=1, 
              help="number of threads used for some bottleneck parts, default=1", metavar="character"),  
  
  make_option(c("-numObsInCluster", "--minimumNumOfElemsInCluster"), type="integer", default=100, 
              help="minimum number of elements in cluster (done for germline), default=100, clustering happens only if number of samples bigger than 3 by number of elements in cluster", metavar="number"),  
  
  make_option(c("-vis", "--visulizationIGV"), action="store_true", default=T, 
              help="if you dont need IGV tracks as output, specify this flag (as printing out IGV tracks slows down the program)"),  
  
  make_option(c("-cloneP", "--clonePenalty"), type="integer", default=200, 
              help="penalty for each additional clone (if you feel that you have some false positive clones, increase this value from default 200)", metavar="number"),  
  
  make_option(c("-d","--debug"), action="store_true", default=FALSE, help="Print debugging information while running.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(paste("We run script located in folder" , opt$folderWithScript, ". All the paths will be calculated realtive to this one. If everything crashes, please, check the correctness of this path first."))




if (is.null(opt$normal) | is.null(opt$bed)) {
  print("You need to specify file with normal coverages and bed file path at least. Here is the help:")
  print_help(opt_parser)
  quit()
}

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
  print("Offtarget files are provided. We try to utilise off-target coverage also (somatic mode).")
  frameworkOff = "offtarget"
}
if (!is.null(opt$normalOfftarget) & !is.null(opt$bedOfftarget) & is.null(opt$tumorOfftarget)) {
  print("Offtarget files are provided. We try to utilise off-target coverage also (germline mode).")
  frameworkOff = "offtargetGermline"
}

frameworkDataTypes = "covdepth"
if (!is.null(opt$bafFolder)) {
  print("Folder with BAFs were provided. Framework swithced to BAF.")
  frameworkDataTypes = "covdepthBAF"
}





no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
cl<-makeCluster(no_cores, type="FORK")
registerDoParallel(cl)



### READING DATA
print(paste("We are started with reading the coverage files and bed files",Sys.time()))
setwd(opt$folderWithScript)
#bedFile <- read.table(opt$bed, stringsAsFactors = F, sep="\t", comment.char="&", header=F)
bedFile <- ReadFileFast(opt$bed, header=F)
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

#normal <- read.table(opt$normal, header=T, stringsAsFactors = F, comment.char="&" )
normal <- ReadFileFast(opt$normal, header=T)
normal = normal
colnames(normal) = cutX(colnames(normal))
if (!startsWith(normal[,1], "chr"))
  normal[,1] <- paste0("chr", normal[,1])
normal <- normal[order(normal[,1], as.numeric(normal[,2])),]
normal <- as.matrix(normal[,opt$colNum:ncol(normal)])
if (length(whichBedIsNA) > 0)
  normal = normal[-whichBedIsNA,]

if (framework == "somatic") {
  #tumor <- read.table(opt$tumor, header=T, stringsAsFactors = F, comment.char="&" )
  tumor <- ReadFileFast(opt$tumor, header=T)
  tumor = tumor
  colnames(tumor) = cutX(colnames(tumor))
  if (!startsWith(tumor[,1], "chr"))
    tumor[,1] <- paste0("chr", tumor[,1])
  tumor <- tumor[order(tumor[,1], as.numeric(tumor[,2])),]
  tumor <- as.matrix(tumor[,opt$colNum:ncol(tumor)])
  if (length(whichBedIsNA) > 0)
    tumor = tumor[-whichBedIsNA,]
}

if (frameworkOff == "offtarget" | frameworkOff == "offtargetGermline") {
  #bedFileOfftarget <- read.table(opt$bedOfftarget, stringsAsFactors = F, sep="\t")
  bedFileOfftarget <- ReadFileFast(opt$bedOfftarget, header=F)
  if (!startsWith(bedFileOfftarget[,1], "chr"))
    bedFileOfftarget[,1] <- paste0("chr", bedFileOfftarget[,1])
  colnames(bedFileOfftarget) <- c("chr.X", "start", "end", "gc")
  bedFileOfftarget <- bedFileOfftarget[order(bedFileOfftarget[,1], as.numeric(bedFileOfftarget[,2])),]
  
  for (i in 1:20) {
    tableOfValues <- table(round(as.numeric(as.character(bedFileOfftarget[,4])) / i, digits = 2) * i)
    if(sum(tableOfValues[which(tableOfValues > 100)]) / sum(tableOfValues) > 0.95) break 
  }
  bedFileOfftarget[,4] <- round(as.numeric(as.character(bedFileOfftarget[,4])) / i, digits = 2) * i
  
  
  #normalOff <- read.table(opt$normalOfftarget, header=T, stringsAsFactors = F, comment.char="&" )
  normalOff <- ReadFileFast(opt$normalOfftarget, header=T)
  colnames(normalOff) = cutX(colnames(normalOff))
  if (!startsWith(normalOff[,1], "chr"))
    normalOff[,1] <- paste0("chr", normalOff[,1])
  normalOff <- normalOff[order(normalOff[,1], as.numeric(normalOff[,2])),]
  normalOff <- as.matrix(normalOff[,opt$colNum:ncol(normalOff)])
  # remain only samples that are in Normal cohort 
  normalOff <- normalOff[,which(colnames(normalOff) %in% colnames(normal))]
  
  #tumorOff <- read.table(opt$tumorOfftarget, header=T, stringsAsFactors = F, comment.char="&" )
  if (frameworkOff == "offtarget") {
    tumorOff <- ReadFileFast(opt$tumorOfftarget, header=T) 
    colnames(tumorOff) = cutX(colnames(tumorOff))
    if (!startsWith(tumorOff[,1], "chr"))
      tumorOff[,1] <- paste0("chr", tumorOff[,1])
    tumorOff <- tumorOff[order(tumorOff[,1], as.numeric(tumorOff[,2])),]
    tumorOff <- as.matrix(tumorOff[,opt$colNum:ncol(tumorOff)])
    # remain only samples that are in Tumor cohort 
    tumorOff <- tumorOff[,which(colnames(tumorOff) %in% colnames(tumor))]
  }
}

print(paste("Started basic quality filtering.",Sys.time()))

rowsToRemove <- cleanDatasetFromLowCoveredFiles(normal, bedFile)
if (length(rowsToRemove) > 0) {
  bedFile <- bedFile[-rowsToRemove,]
  normal <- normal[-rowsToRemove,]
  if (framework == "somatic")
    tumor <- tumor[-rowsToRemove,]
}
if (frameworkOff == "offtarget" | frameworkOff == "offtargetGermline") {
  rowsToRemove <- cleanDatasetFromLowCoveredFiles(normalOff, bedFile)
  if (length(rowsToRemove) > 0) {
    normalOff = normalOff[-rowsToRemove,]
    if (frameworkOff == "offtarget")
      tumorOff = tumorOff[-rowsToRemove,]
    bedFileOfftarget = bedFileOfftarget[-rowsToRemove,]
  }
}

### GC CONTENT NORMALIZATION

if (framework == "somatic") {
  pairs <- read.table(opt$pair, sep=",", stringsAsFactors = F)
  pairs <- data.frame(pairs, ncol=2)
  pairs <- unique(pairs)
}


## CHECK INPUT VALIDITY
if (!is.null(opt$normalSample)) {
  if (!opt$normalSample %in% colnames(normal)) {
    print(paste("Sorry! We have not found sample with name", opt$normalSample, "in our normal cohort! Please check that sample name matches. Header of cohort is:"))
    print(colnames(normal))
  }
  stopifnot(opt$normalSample %in% colnames(normal))
  print(paste("We have found sample with name", opt$normalSample, "in our normal.cov file! We will try to analyse it.", 
              "Remember - we still have to calculate all the parameters (the most time consuming step)..."))
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
  print(paste("We are reading BAF files. It may take time - especially if you have a lot of SNV positions.", Sys.time()))
  setwd(opt$folderWithScript)

    source(file.path(opt$folderWithScript, "somatic", "bafSegmentation.R"), local=T)

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
  if (length(allowedChromsBaf) == 0) {
    print("Apparently none of your baf files match with sample pairs you've provided. We can not use any bafs from now on and rely only on coverage.")
    frameworkDataTypes = "covdepth"
    rm(allowedChromsBaf)
    rm(bAlleleFreqsAllSamples)
  }
}
setwd(opt$folderWithScript)



### ON TARGET GC NORMALIZATION
print(paste("Normalization with GC and length starts.", Sys.time()))
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
if (frameworkOff == "offtarget" | frameworkOff == "offtargetGermline") {
  lst <- gc_and_sample_size_normalise(bedFileOfftarget, normalOff)
  normalOff <- lst[[1]]
  if (frameworkOff == "offtarget") {
    if (frameworkDataTypes == "covdepthBAF") {
      lst <- gc_and_sample_size_normalise(bedFileOfftarget, tumorOff, allowedChroms=allowedChromsBaf)
    } else {
      lst <- gc_and_sample_size_normalise(bedFileOfftarget, tumorOff)
    }
  }
  if (frameworkOff == "offtarget") {
    tumorOff <- lst[[1]]
  }
  bedFileOfftarget <- lst[[2]]
}

### BED FILE OFFTARGET MAY NOT CONTAIN COLUMN WITH GENES
if (ncol(bedFile)  == 4) {
  bedFile <- cbind(bedFile, rep(0, nrow(bedFile)))
  colnames(bedFile) <- colnames(bedFile)
}
if (frameworkOff == "offtarget" | frameworkOff == "offtargetGermline") {
  if (ncol(bedFileOfftarget)  == 4) {
    bedFileOfftarget <- cbind(bedFileOfftarget, rep(0, nrow(bedFileOfftarget)))
    colnames(bedFileOfftarget) <- colnames(bedFile)
  }
}

# FILTER LOW COVERED REGIONS
if (framework == "somatic") {
  regionsToFilerOutOn = findRegionsToFilerOutDueSystematicallyLowCoverage(normal, tumor)
} else {
  regionsToFilerOutOn = findRegionsToFilerOutDueSystematicallyLowCoverage(normal)
}

if (length(regionsToFilerOutOn)>0) {
	normal = normal[-regionsToFilerOutOn,] + 10**-20
	if (framework == "somatic") {
	  tumor = tumor[-regionsToFilerOutOn,] + 10**-20
	}
	bedFile = bedFile[-regionsToFilerOutOn,]
}


if (frameworkOff == "offtarget" | frameworkOff == "offtargetGermline") {
  if (frameworkOff == "offtarget") {
    regionsToFilerOutOff = findRegionsToFilerOutDueSystematicallyLowCoverage(normal, tumor)
  } else {
    regionsToFilerOutOff = findRegionsToFilerOutDueSystematicallyLowCoverage(normal)
  }
  if (length(regionsToFilerOutOff) > 0){
    normalOff = normalOff[-regionsToFilerOutOff,] + 10**-20
    if (frameworkOff == "offtarget")
      tumorOff = tumorOff[-regionsToFilerOutOff,] + 10**-20
    bedFileOfftarget = bedFileOfftarget[-regionsToFilerOutOff,]
  }
}

### EXTRACTING INFORMATION FROM BED
bordersOfChroms <- getBordersOfChromosomes(bedFile)




stopCluster(cl)
### PROCESSING OF GERMLINE VARIANTS


setwd(opt$folderWithScript)
source("./germline/helpersGermline.R")

frameworkTrios = "single"
if (!is.null(opt$triosFile)) {
  frameworkTrios = "trios"
}

if (frameworkTrios == "trios") {
  trios <- read.table(opt$triosFile, sep=",", stringsAsFactors = F)
  trios <- data.frame(trios)
  trios <- unique(trios)
  colnames(trios) <- c("Kid","Mother","Father")
}

print(paste("We start to cluster your data (you will find a plot if clustering is possible in your output directory)", opt$out, Sys.time()))
clustering <- returnClustering(as.numeric(opt$minimumNumOfElemsInCluster))

print(paste("Processing of germline variants started (we need to do it as an additional step for Saomtic calling since we need to know at least genders).", Sys.time()))

orderOfBedFile <- order(bedFile[,1], as.numeric(bedFile[,2]))
bedFile = bedFile[orderOfBedFile,]
normal = normal[orderOfBedFile,]
coverageAllSamples <- sqrt(as.matrix(normal))

if (frameworkOff == "offtargetGermline") {
  orderOfBedFileOff <- order(bedFileOfftarget[,1], as.numeric(bedFileOfftarget[,2]))
  bedFileOfftarget = bedFileOfftarget[orderOfBedFileOff,]
  normalOff = normalOff[orderOfBedFileOff,]
  coverageAllSamplesOff <- sqrt(as.matrix(normalOff))
}

print(paste("Gender estimation started", Sys.time()))
genderOfSamplesCohort <- Determine.gender(sqrt(coverageAllSamples), bedFile)
print(genderOfSamplesCohort)
print(paste("Gender succesfully determined. Plot is written in your results directory:", opt$out, Sys.time()))
  
if (framework == "germline") {
  for (cluster in unique(clustering)) {
    if (cluster == -1) {
      print(paste("Samples from trio mode that are presented in trios.txt but do not have a full family in file", opt$normal , "will be excluded."))
      print(colnames(normal)[which(clustering == -1)])
      next
    }
    
    # We create cluster for parallel computation each time we run germline analysis
    no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
    cl<-makeCluster(no_cores, type="FORK")
    registerDoParallel(cl)
    
    
    samplesToAnalyse = which(clustering == cluster)
    coverage <- coverageAllSamples[,samplesToAnalyse]
    genderOfSamples = genderOfSamplesCohort[samplesToAnalyse]
    if (frameworkOff == "offtargetGermline") {
      samplesToAnalyseOff = which(colnames(coverageAllSamplesOff) %in% colnames(coverage))
      coverageOff <- coverageAllSamplesOff[,samplesToAnalyseOff]
      genderOfSamplesOff <- sapply(1:ncol(coverageOff), function(i) {return(genderOfSamplesCohort[which(colnames(normal) == colnames(coverageOff)[i])])})
    }
    
    
    
    #clusterExport(cl, c('EstimateModeSimple', 'bedFile', 'genderOfSamples', "lehmanHodges", 'Qn'))
    
    
    
    
    
    autosomes <- which(!bedFile[,1] %in% c("chrX", "chrY", "X", "Y"))
    sdsOfGermlineSamples <- apply(coverage[autosomes,], 2, determineSDsOfGermlineSample)
    
    #medians <- parSapply(cl=cl, 1:nrow(coverage), function(i) {EstimateModeSimple(coverage[i,], bedFile[i,1], FindRobustMeanAndStandardDeviation)})
    print(paste("We start estimation of parameters of germline cohort. It may take some time. Cluster of samples being analysed:", cluster, Sys.time()))
    mediansAndSds = calculateLocationAndScale(bedFile, coverage, genderOfSamples, sdsOfGermlineSamples, autosomes)
    medians = as.numeric(mediansAndSds[,1])
    sdsOfProbes = trimValues(as.numeric(mediansAndSds[,2]), 0.01)
    coverage.normalised = sweep(coverage, 1, medians, FUN="/")
    
    if (frameworkOff == "offtargetGermline") {
      autosomesOff <- which(!bedFileOfftarget[,1] %in% c("chrX", "chrY", "X", "Y"))
      sdsOfGermlineSamplesOff <- apply(coverageOff[autosomesOff,], 2, determineSDsOfGermlineSample)
      mediansAndSdsOff = calculateLocationAndScale(bedFileOfftarget, coverageOff, genderOfSamplesOff, sdsOfGermlineSamplesOff, autosomesOff)
      mediansOff = as.numeric(mediansAndSdsOff[,1])
      sdsOfProbesOff = trimValues(as.numeric(mediansAndSdsOff[,2]), 0.01)
      coverage.normalised.off = sweep(coverageOff, 1, mediansOff, FUN="/")
    }
    
    
    stopCluster(cl)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if (!is.null(opt$triosFile)) {
      source("./trios/germlineTrioSolver.R",local=TRUE)
    } else {
      source("./germline/germlineSolver.R",local=TRUE)
    }
  }
  if (framework == "germline" | !is.null(opt$triosFile)) quit()
}




print(paste("Processing of somatic variants started.", Sys.time()))

no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
cl<-makeCluster(no_cores, type="FORK")
registerDoParallel(cl)



setwd(opt$folderWithScript)
source("./somatic/somaticSolver.R",local=TRUE)
stopCluster(cl)
