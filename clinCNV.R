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
  
  make_option(c("-sg", "--scoreG"), type="double", default="20", 
              help="minimum threshold for significance germline variants", metavar="number"),
  
  make_option(c("-lg", "--lengthG"), type="integer", default="2", 
              help="minimum threshold for length of germline variants", metavar="number"),
  
  make_option(c("-ss", "--scoreS"), type="double", default="100", 
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
              help="penalty for each additional clone (if you feel that you have some false positive clones, increase this value from default 200)"),  
  
  make_option(c("-purityS", "--purityStep"), type="double", default=2.5, 
              help="step of purity we investigate (from 5% to 100% with the step you specify, default=2.5)", metavar="number"),  
  
  make_option(c("-dfStudent", "--degreesOfFreedomStudent"), type="integer", default=1000, 
              help="number of degrees of freedom of Student's distribution for somatic analysis (a lot of outliers => reduce the default value of 1000 to e.g. 10)"),  
  
  make_option(c("-polymC", "--polymorphicCalling"), type="character", default="NO", 
              help="should calling of polymorphic regions be performed, YES = calling is performed, NO = no polymorphic calling (default), any other string = mCNVs taken from the file with that path (it must have at least 3 columns chrom-start-end)"),  
  
  make_option(c("-mosaic", "--mosaicism"), action="store_true", default=F, 
              help="if mosaic calling should be performed"),  
  
  make_option(c("-d","--debug"), action="store_true", default=FALSE, help="Print debugging information while running.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt$folderWithScript = normalizePath(opt$folderWithScript)
print(paste("We run script located in folder" , opt$folderWithScript, ". All the paths will be calculated realtive to this one. If everything crashes, please, check the correctness of this path first."))

### TESTING PART
opt$bed = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/ssSC_v4.annotated.bed"
opt$tumor = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/tumor_ontarget_v4.cov"
opt$normal = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/ontarget_v4.cov"
opt$colNum = 4
opt$pair = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/pairsNew.txt"
opt$out = "/Users/gdemidov/Tuebingen/clinCNV_dev/results"
opt$reanalyseCohort = F
opt$bedOfftarget = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/annotated_offtarget_v4.bed"
opt$tumorOfftarget = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/tumor_offtarget_v4.cov"
opt$normalOfftarget = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/offtarget_v4.cov"
opt$bafFolder = "/Users/gdemidov/Tuebingen/somatic_CNVs/Somatic/baf"
opt$folderWithScript = "/Users/gdemidov/Tuebingen/clinCNV_dev_new/ClinCNV/"

if (is.null(opt$normal) | is.null(opt$bed)) {
  print("You need to specify file with normal coverages and bed file path at least. Here is the help:")
  print_help(opt_parser)
  quit()
}

setwd(opt$folderWithScript)
source(paste0(opt$folderWithScript, "/generalHelpers.R"))

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

if (opt$polymorphicCalling == "YES") {
  print("You've choosen to detect polymorphic regions with the help of our tool - great choice!")
}

if (opt$mosaicism) {
  print("You suspect your samples to be mosaic - hmmm, we will check this out...(but the mosaic CN change should not be > 1 copy different from default")
}



no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
cl<-makeCluster(no_cores, type="FORK")
registerDoParallel(cl)



### READING DATA
print(paste("We are started with reading the coverage files and bed files",Sys.time()))
setwd(opt$folderWithScript)
bedFile <- ReadFileFast(opt$bed, header=F)
if (!startsWith(bedFile[,1], "chr"))
  bedFile[,1] <- paste0("chr", bedFile[,1])
if (ncol(bedFile)  == 4) {
  bedFile <- cbind(bedFile, rep(0, nrow(bedFile)))
  colnames(bedFile) <- colnames(bedFile)
}
colnames(bedFile) <- c("chr.X", "start", "end", "gc", "genes")
bedFile <- bedFile[order(bedFile$chr.X, as.numeric(bedFile$start)),]
presentedChromsOn = unique(bedFile[,1])
numberOfElemsInEachChromosome = sapply(1:length(presentedChromsOn), function(i) {
   if (length(which(bedFile[,1] == presentedChromsOn[i])) > 10) {
     return(T)
   } else {
     return(F)
   }
})



for (i in 1:20) {
  tableOfValues <- table(round(as.numeric(as.character(bedFile[,4])) / i, digits = 2) * i)
  if(sum(tableOfValues[which(tableOfValues > 100)]) / sum(tableOfValues) > 0.95) break 
}
bedFile[,4] <- round(as.numeric(as.character(bedFile[,4])) / i, digits = 2) * i
whichBedIsNA <- which(is.na(bedFile[,4]) | bedFile[,3] - bedFile[,2] < 80 | which(!bedFile[,1] %in% presentedChromsOn[numberOfElemsInEachChromosome]))
if (length(whichBedIsNA) > 0)
  bedFile = bedFile[-whichBedIsNA,]

normal <- ReadFileFast(opt$normal, header=T)
#colnames(normal) = c("chr","start","end", 1:(ncol(normal) - 3))
colnames(normal) = cutX(colnames(normal))
if (!startsWith(normal[,1], "chr"))
  normal[,1] <- paste0("chr", normal[,1])
normal <- normal[order(normal[,1], as.numeric(normal[,2])),]
normal <- as.matrix(normal[,opt$colNum:ncol(normal)])
normal = checkForDuplicatesAndRemove(normal, opt$normalSample)
if (length(whichBedIsNA) > 0)
  normal = normal[-whichBedIsNA,]

numberOfRowsBeforeAllTheFiltrationNormal = nrow(normal)

if (framework == "somatic") {
  tumor <- ReadFileFast(opt$tumor, header=T)
  colnames(tumor) = cutX(colnames(tumor))
  if (!startsWith(tumor[,1], "chr"))
    tumor[,1] <- paste0("chr", tumor[,1])
  tumor <- tumor[order(tumor[,1], as.numeric(tumor[,2])),]
  tumor <- as.matrix(tumor[,opt$colNum:ncol(tumor)])
  tumor = checkForDuplicatesAndRemove(tumor, opt$tumorSample)
  if (length(whichBedIsNA) > 0)
    tumor = tumor[-whichBedIsNA,]
}


if (frameworkOff == "offtarget" | frameworkOff == "offtargetGermline") {
  bedFileOfftarget <- ReadFileFast(opt$bedOfftarget, header=F)
  if (!startsWith(bedFileOfftarget[,1], "chr"))
    bedFileOfftarget[,1] <- paste0("chr", bedFileOfftarget[,1])
  ### BED FILE OFFTARGET MAY NOT CONTAIN COLUMN WITH GENES
  
  if (frameworkOff == "offtarget" | frameworkOff == "offtargetGermline") {
    if (ncol(bedFileOfftarget)  == 4) {
      bedFileOfftarget <- cbind(bedFileOfftarget, rep(0, nrow(bedFileOfftarget)))
      colnames(bedFileOfftarget) <- colnames(bedFile)
    }
  }
  colnames(bedFileOfftarget) <- c("chr.X", "start", "end", "gc", "genes")
  bedFileOfftarget <- bedFileOfftarget[order(bedFileOfftarget[,1], as.numeric(bedFileOfftarget[,2])),]
  
  for (i in 1:20) {
    tableOfValues <- table(round(as.numeric(as.character(bedFileOfftarget[,4])) / i, digits = 2) * i)
    if(sum(tableOfValues[which(tableOfValues > 100)]) / sum(tableOfValues) > 0.95) break 
  }
  bedFileOfftarget[,4] <- round(as.numeric(as.character(bedFileOfftarget[,4])) / i, digits = 2) * i
  presentedChromsOff = unique(bedFileOfftarget[,1])
  numberOfElemsInEachChromosomeOff = sapply(1:length(presentedChromsOff), function(i) {
    if (length(which(bedFileOfftarget[,1] == presentedChromsOff[i])) > 10) {
      return(T)
    } else {
      return(F)
    }
  })
  whichBedOffIsNA <- which(is.na(bedFileOfftarget[,4]) | (!bedFileOfftarget[,1] %in% presentedChromsOff[numberOfElemsInEachChromosomeOff]))
  bedFileOfftarget = bedFileOfftarget[-whichBedOffIsNA,]
  
  normalOff <- ReadFileFast(opt$normalOfftarget, header=T)
  colnames(normalOff) = cutX(colnames(normalOff))
  if (!startsWith(normalOff[,1], "chr"))
    normalOff[,1] <- paste0("chr", normalOff[,1])
  
  normalOff <- normalOff[order(normalOff[,1], as.numeric(normalOff[,2])),]
  normalOff <- as.matrix(normalOff[,opt$colNum:ncol(normalOff)])
  # remain only samples that are in Normal cohort 
  normalOff <- normalOff[,which(colnames(normalOff) %in% colnames(normal))]
  normalOff <- checkForDuplicatesAndRemove(normalOff, opt$normalSample)
  normalOff <- normalOff[-whichBedOffIsNA,]
  
  if (frameworkOff == "offtarget") {
    tumorOff <- ReadFileFast(opt$tumorOfftarget, header=T) 
    colnames(tumorOff) = cutX(colnames(tumorOff))
    if (!startsWith(tumorOff[,1], "chr"))
      tumorOff[,1] <- paste0("chr", tumorOff[,1])
    tumorOff <- tumorOff[order(tumorOff[,1], as.numeric(tumorOff[,2])),]
    tumorOff <- as.matrix(tumorOff[,opt$colNum:ncol(tumorOff)])
    # remain only samples that are in Tumor cohort 
    tumorOff <- tumorOff[,which(colnames(tumorOff) %in% colnames(tumor))]
    tumorOff <- checkForDuplicatesAndRemove(tumorOff, opt$tumorSample)
    tumorOff <- tumorOff[-whichBedOffIsNA,]
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
print(paste("Amount of regions after filtering of 0-covered regions", round(100 * nrow(normal) / numberOfRowsBeforeAllTheFiltrationNormal, digits = 3) ) )
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
if (max(bedFile[,3] - bedFile[,2]) / min(bedFile[,3] - bedFile[,2]) > 16) {
  lengthBasedNorm = T
  normal <- lengthBasedNormalization(normal, bedFile)
} else {
  lengthBasedNorm = F
}
  
lst <- gc_and_sample_size_normalise(bedFile, normal)
normal <- lst[[1]]
if (framework == "somatic") {
  if (frameworkDataTypes == "covdepthBAF") {
    if (lengthBasedNorm)
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
    tumorOff <- lst[[1]]
  }
  bedFileOfftarget <- lst[[2]]
}

print(paste("Amount of regions after GC-extreme filtering", round(100 * nrow(normal) / numberOfRowsBeforeAllTheFiltrationNormal, digits = 3) ) )



# FILTER LOW COVERED REGIONS
if (framework == "somatic") {
  regionsToFilerOutOn = findRegionsToFilerOutDueSystematicallyLowCoverage(normal, tumor)
} else {
  regionsToFilerOutOn = findRegionsToFilerOutDueSystematicallyLowCoverage(normal)
}
print(paste("Amount of regions after Systematically Low Covered regions filtering", round(100 * nrow(normal) / numberOfRowsBeforeAllTheFiltrationNormal, digits = 3) ) )

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
source(paste0(opt$folderWithScript, "/germline/helpersGermline.R"))

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
clusteringList <- returnClustering(as.numeric(opt$minimumNumOfElemsInCluster))
clustering = clusteringList[[1]]
outliersByClusteringCohort = clusteringList[[2]]

  
orderOfBedFile <- order(bedFile[,1], as.numeric(bedFile[,2]))
bedFile = bedFile[orderOfBedFile,]
normal = normal[orderOfBedFile,]
coverageAllSamples <- sqrt(as.matrix(normal))



print(paste("Gender estimation started", Sys.time()))
genderOfSamplesCohort <- Determine.gender(coverageAllSamples, bedFile)
print(genderOfSamplesCohort)
print(paste("Gender succesfully determined. Plot is written in your results directory:", opt$out, Sys.time()))

if (framework == "germline") {
  gc()
  print(paste("Processing of germline variants started (we need to do it as an additional step for Saomtic calling since we need to know at least genders).", Sys.time()))
  

  
  if (frameworkOff == "offtargetGermline") {
    orderOfBedFileOff <- order(bedFileOfftarget[,1], as.numeric(bedFileOfftarget[,2]))
    bedFileOfftarget = bedFileOfftarget[orderOfBedFileOff,]
    normalOff = normalOff[orderOfBedFileOff,]
    coverageAllSamplesOff <- sqrt(as.matrix(normalOff))
  }
  
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
    if (!is.null(opt$normalSample)) {
      if (!opt$normalSample %in% colnames(coverageAllSamples)[samplesToAnalyse]) {
        next
      }
    }
    
    coverage <- coverageAllSamples[,samplesToAnalyse]
    genderOfSamples = genderOfSamplesCohort[samplesToAnalyse]
    outliersByClustering = outliersByClusteringCohort[samplesToAnalyse]
    if (frameworkOff == "offtargetGermline") {
      samplesToAnalyseOff = which(colnames(coverageAllSamplesOff) %in% colnames(coverage))
      coverageOff <- coverageAllSamplesOff[,samplesToAnalyseOff]
      genderOfSamplesOff <- sapply(1:ncol(coverageOff), function(i) {return(genderOfSamplesCohort[which(colnames(normal) == colnames(coverageOff)[i])])})
    }
    
    
    
    #clusterExport(cl, c('EstimateModeSimple', 'bedFile', 'genderOfSamples', "lehmanHodges", 'Qn'))
    
    
    polymorphicRegions = NULL
    ### HERE THE POLYMORPHIC REGIONS DETECTION GOES
    if (opt$polymorphicCalling == "YES") {
      source(paste0(opt$folderWithScript, "/germline/mCNVsDetection.R"),local=TRUE)
      if (nrow(copyNumberForReporting) > 0)
      polymorphicRegions = copyNumberForReporting
    }
    if (opt$polymorphicCalling != "YES" & opt$polymorphicCalling != "NO") {
      print("Since polymorphicCalling option was not YES or NO, we interpret the value as a path to file")
      print(opt$polymorphicCalling)
      print("If ClinCNV fails now, then check the path to the file! It should have at least 3 columns: chrom, start, end.")
      polymorphicRegions = read.table(opt$polymorphicCalling, header=T)
    }
    
    
    autosomes <- which(!bedFile[,1] %in% c("chrX", "chrY", "X", "Y"))

    #medians <- parSapply(cl=cl, 1:nrow(coverage), function(i) {EstimateModeSimple(coverage[i,], bedFile[i,1], FindRobustMeanAndStandardDeviation)})
    print(paste("We start estimation of parameters of germline cohort. It may take some time. Cluster of samples being analysed:", cluster, Sys.time()))
    mediansAndSds = calculateLocationAndScale(bedFile, coverage, genderOfSamples, autosomes)
    coverage.normalised = mediansAndSds[[1]]
    
    

    
    sdsOfProbes = trimValues(as.numeric(mediansAndSds[[2]][,2]), 0.01)
    sdsOfGermlineSamples = mediansAndSds[[3]]
    filterOutRegionsWithSmallMedians <- which(mediansAndSds[[2]][,1] > 0.3)
    sdsOfProbes = sdsOfProbes[filterOutRegionsWithSmallMedians]
    bedFileFiltered = bedFile[filterOutRegionsWithSmallMedians,]
    coverage.normalised = coverage.normalised[filterOutRegionsWithSmallMedians,]
    print(paste("Amount of regions after low covered regions filtering (median < 0.3) in cluster", cluster, ":", round(100 * nrow(coverage.normalised) / numberOfRowsBeforeAllTheFiltrationNormal, digits = 3) ) )
    
    if (frameworkOff == "offtargetGermline") {
      
      autosomesOff = which(!bedFileOfftarget[,1] %in% c("chrX", "chrY"))
      mediansAndSdsOff = calculateLocationAndScale(bedFileOfftarget, coverageOff, genderOfSamplesOff, autosomesOff)
      coverage.normalised.off = mediansAndSdsOff[[1]]
      sdsOfProbesOff = trimValues(as.numeric(mediansAndSdsOff[[2]][,2]), 0.01)
      sdsOfGermlineSamplesOff = mediansAndSdsOff[[3]]
      filterOutRegionsWithSmallMedians <- which(mediansAndSdsOff[[2]][,1] > 0.3)
      sdsOfProbesOff = sdsOfProbesOff[filterOutRegionsWithSmallMedians]
      bedFileFilteredOfftarget = bedFileOfftarget[filterOutRegionsWithSmallMedians,]
      coverage.normalised.off = coverage.normalised.off[filterOutRegionsWithSmallMedians,]
    }
    
    
    stopCluster(cl)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if (!is.null(opt$triosFile)) {
      source(paste0(opt$folderWithScript,"/trios/germlineTrioSolver.R"),local=TRUE)
    } else {
      source(paste0(opt$folderWithScript, "/germline/germlineSolver.R"),local=TRUE)
    }
  }
  if (framework == "germline" | !is.null(opt$triosFile)) quit()
}




print(paste("Processing of somatic variants started.", Sys.time()))



setwd(opt$folderWithScript)

for (cluster in unique(clustering)) {
    print(paste("Working on somatic samples, cluster", cluster, Sys.time()))
  
    samplesToAnalyse = which(clustering == cluster)
    genderOfSamples = genderOfSamplesCohort[samplesToAnalyse]
    tmpNormal = normal[,which(clustering == cluster)]
    if (!is.null(opt$normalSample) & !is.null(opt$tumorSample)) {
      if (!opt$normalSample %in% colnames(tmpNormal)) {
        next
      }
    }
    source(paste0(opt$folderWithScript, "/somatic/somaticSolver.R"),local=TRUE)
}


