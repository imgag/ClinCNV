#!/usr/bin/env Rscript
set.seed(100)
options(warn=-1)
clincnvVersion = paste0("ClinCNV version: v1.18.1")

## CHECK R VERSION
if (!( (as.numeric(version$major) >= 3 & as.numeric(version$minor) > 2.0) |  as.numeric(version$major) >= 4) ) {
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
library(mclust)
library(R.utils)
library(umap)
library(dbscan)


Rcpp_global = "Rcpp" %in% rownames(installed.packages())
if (Rcpp_global) {library("Rcpp")}


initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)


## DETERMINE THE PATH TO THE SCRIPT AUTOMATICALLY
current_working_dir <- script.basename


option_list = list(
  make_option("--normal", type="character", default=NULL, 
              help="path to table with normal coverages"),
  
  make_option("--tumor", type="character", default=NULL, 
              help="path to table with tumor coverages"),
  
  make_option("--normalOfftarget", type="character", default=NULL, 
              help="path to table with normal offtarget coverages"),
  
  make_option("--tumorOfftarget", type="character", default=NULL, 
              help="path to table with tumor offtarget coverages"),
  
  make_option("--out", type="character", default="./result/", 
              help="output folder path [default= %default]"),
  
  make_option("--pair", type="character", default="pairs.txt", 
              help="file with pairing information, 1st column = tumor, 2nd column = normal [default= %default]"),
  
  make_option("--bed", type="character", default=NULL, 
              help="bed file with panel description (chr \t start \t end \t gc_content \t annotation). has to use same notation as .cov files."),
  
  make_option("--bedOfftarget", type="character", default=NULL, 
              help="offtarget bed file with panel description (chr \t start \t end \t gc_content \t annotation). has to use same notation as .cov files."),
  
  make_option("--colNum", type="integer", default=4, 
              help="column where coverages start"),
  
  make_option("--folderWithScript", type="character", default=current_working_dir, 
              help="folder where you put script"),
  
  make_option("--reanalyseCohort", action="store_false", 
              help="if specified, reanalyses whole cohort [default= %default]"),
  
  make_option("--scoreG", type="double", default="20", 
              help="minimum threshold for significance germline variants"),
  
  make_option("--lengthG", type="integer", default="2", 
              help="minimum threshold for length of germline variants"),
  
  make_option("--scoreS", type="double", default="100", 
              help="minimum threshold for significance somatic variants"),
  
  make_option("--lengthS", type="integer", default="9", 
              help="minimum threshold for length of somatic variants"),
  
  make_option("--maxNumGermCNVs", type="integer", default="10000", 
              help="maximum number of germline CNVs allowed (increase thresholds if does not meet criteria)"),
  
  make_option("--maxNumSomCNAs", type="integer", default="10000", 
              help="maximum number of somatic CNAs allowed (increase thresholds if does not meet criteria)"),
  
  make_option("--maxNumIter", type="integer", default=3, 
              help="maximum number of iterations of variant calling"),
  
  make_option("--bafFolder", type="character", default=NULL, 
              help="folder where you put BAF frequencies (one per normal, one per tumor sample)"),
  
  make_option("--normalSample", type="character", default=NULL, 
              help="name of normal sample to analyse (if only one sample has to be analysed)"),
  
  make_option("--tumorSample", type="character", default=NULL, 
              help="name of tumor sample to analyse (if only one sample has to be analysed, normal has to be provided too)"),
  
  make_option("--triosFile", type="character", default=NULL, 
              help="file with information about trios, child-father-mother"),
  
  make_option("--fdrGermline", type="integer", default=0, 
              help="number of iterations for FDR check (more - better, but slower, 0 = no FDR correction)"),
  
  make_option("--numberOfThreads", type="integer", default=1, 
              help="number of threads used for some bottleneck parts, default=1"),  
  
  make_option("--minimumNumOfElemsInCluster", type="integer", default=10000, 
              help="minimum number of elements in cluster (done for germline), default=100, clustering happens only if number of samples bigger than 3 by number of elements in cluster", metavar="number"),  
  
  make_option("--visulizationIGV", action="store_true", default=T, 
              help="if you dont need IGV tracks as output, specify this flag (as printing out IGV tracks slows down the program)"),  
  
  make_option("--clonePenalty", type="integer", default=300, 
              help="penalty for each additional clone (if you feel that you have some false positive clones, increase this value from default 300)"),  
  
  make_option("--purityStep", type="double", default=2.5, 
              help="step of purity we investigate (from 5% to 100% with the step you specify, default=2.5)", metavar="number"),  
  
  make_option("--degreesOfFreedomStudent", type="integer", default=1000, 
              help="number of degrees of freedom of Student's distribution for somatic analysis (a lot of outliers => reduce the default value of 1000 to e.g. 10)"),  
  
  make_option("--polymorphicCalling", type="character", default="NO", 
              help="should calling of polymorphic regions be performed, YES = calling is performed, NO = no polymorphic calling (default), any other string = mCNVs taken from the file with that path (it must have at least 3 columns chrom-start-end)"),  
  
  make_option("--mosaicism", action="store_true", default=F, 
              help="if mosaic calling should be performed"),  
  
  make_option("--minimumPurity", type="double", default=5, 
              help="minimum purity for somatic samples"),  
  
  make_option("--superRecall", type="double", default=10000, 
              help="Super recall mode - after calling normal CNVs it tries to find CNVs with any length that are better than pre-specified threshold"),  
  
  make_option("--clonalityForChecking", type="double", default=0.4, 
              help="Starting from which clonality BAF-based QC-control has to be applied (no allelic balanced variants with smaller purity will be detected!)"),  
  
  make_option("--shiftToTry", type="integer", default=1, 
              help="change only if you have a sample with lots of allelic imbalance (if you think that the diploid baseline should be different, number of options for choosing will be provided during calling)"),  
  
  make_option("--filterStep", type="integer", default=1, 
              help="This value indicates if ClinCNV should perform QC internally (starting from threshold specified by --clonalityForChecking). Value 0 means no, value 1 - only for finding clonality, value 2 - for clonality and final calls too"),  
  
  make_option("--guideBaseline", type="character", default=NULL, 
              help="For complex samples with potential whole-genome duplication - string denoting which region you suspect to be diploid so tool will take it is a baseline (format chrN:12345-67890)"),  
  
  make_option("--notComplexTumor", action="store_true", default=F, 
              help="Sometimes some CNAs happen in the same region twice and leave the signature unrecognizable by simple models. Specify this flag if you don't want the 2nd CNAs to be recognized by ClinCNV."),  
  
  make_option("--pnealtyHigherCopy", type="double", default=1, 
              help="How big should be penalty for higher copy? This is penalty for each additional copy, one per CNV. (smaller values: more big copy-number allowed, lower clonal cancer cell fraction)"),  
  
  make_option("--pnealtyHigherCopyOneSegment", type="double", default=0.01, 
              help="How big should be penalty for higher copy? This is penalty for each additional copy, one per region in CNV. (smaller values: more big copy-number allowed, lower clonal cancer cell fraction)"),  
  
  make_option("--par", type="character", default="NO", 
              help="coordinates of chrX paralogous regions (format chrX:60001-2699520;chrX:154931044-155260560 )"),  
  
  make_option("--sex", type="character", default="", 
              help="override the sample's gender (active only when you specify --normalSample flag)"),  
  
  make_option("--clusterProvided", type="character", default=NULL, 
              help="Use external clustering (file with lines, sample ID \t cluster ID"),
  
  make_option("--maximumSomaticCN", type="integer", default=30, 
              help="The highest allowed somatic copy-number (higher = more accurate, but slower), [default= %default]"),
  
  make_option("--onlyTumor", action="store_true", default=F, 
              help="if tumor only calling is to be performed"),  
  
  make_option("--noPlot", action="store_true", default=F, 
              help="Do not perform additional plotting"), 
  
  make_option("--hg38", action="store_true", default=F, 
              help="Work with hg38 cytobands switch"),  
  
  make_option(c("-d","--debug"), action="store_true", default=FALSE, help="Print debugging information while running.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt$folderWithScript = normalizePath(opt$folderWithScript)
print(paste("We run script located in folder" , opt$folderWithScript, ". Please, specify ABSOLUTE paths, relative paths do not work for every machine. If everything crashes, please, check the correctness of this path first."))





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

if (opt$onlyTumor) {
  print("You want to detect Tumor Only CN changes - please not that only simple tumors can be analyzed with this flag...")
}




#no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
#cl<-makeCluster(no_cores, type="FORK")
#registerDoParallel(cl)


cl = NULL

print("START cluster allocation.")
no_cores <- min(detectCores() - 1, as.numeric(opt$numberOfThreads))
cl = makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
print("Cluster allocated.")
print("END cluster allocation.")

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
if (!opt$colNum == 1)
  bedFile <- bedFile[order(bedFile$chr.X, as.numeric(bedFile$start)),]
presentedChromsOn = unique(bedFile[,1])
numberOfElemsInEachChromosome = sapply(1:length(presentedChromsOn), function(i) {
  if (length(which(bedFile[,1] == presentedChromsOn[i])) > 3) {
    return(T)
  } else {
    return(F)
  }
})

bedPositionsThatWillBeFiltered = matrix(nrow=0, ncol=ncol(bedFile) + 1)

for (i in 1:20) {
  tableOfValues <- table(round(as.numeric(as.character(bedFile[,4])) / i, digits = 2) * i)
  if(sum(tableOfValues[which(tableOfValues > 25)]) / sum(tableOfValues) > 0.95) break 
}
bedFile[,4] <- round(as.numeric(as.character(bedFile[,4])) / i, digits = 2) * i
whichBedIsNA <- which(is.na(bedFile[,4]) | bedFile[,3] - bedFile[,2] < 50 | (!bedFile[,1] %in% presentedChromsOn[numberOfElemsInEachChromosome]))
bedPositionsThatWillBeFiltered = rbind(bedPositionsThatWillBeFiltered, cbind(bedFile[whichBedIsNA,], rep("TooShortOrNA", length(whichBedIsNA))))
colnames(bedPositionsThatWillBeFiltered)[ncol(bedPositionsThatWillBeFiltered)] = "Description"
if (length(whichBedIsNA) > 0)
  bedFile = bedFile[-whichBedIsNA,]

normal <- ReadFileFast(opt$normal, header=T)
#colnames(normal) = c("chr","start","end", 1:(ncol(normal) - 3))
colnames(normal) = cutX(colnames(normal))
if (!startsWith(normal[,1], "chr") & !opt$colNum == 1)
  normal[,1] <- paste0("chr", normal[,1])
if (!opt$colNum == 1)
  normal <- normal[order(normal[,1], as.numeric(normal[,2])),]

if (length(whichBedIsNA) > 0)
  normal = normal[-whichBedIsNA,]

if (!opt$colNum == 1) {
  if (!checkBedAndCoverageValidity(bedFile, normal)) {
    quit()
  }
}

normal <- as.matrix(normal[,opt$colNum:ncol(normal)])
normal = checkForDuplicatesAndRemove(normal, opt$normalSample)

avgDepthNormalOn = determineAverageDepth(normal, bedFile)
numberOfRowsBeforeAllTheFiltrationNormal = nrow(normal)




if (framework == "somatic") {
  tumor <- ReadFileFast(opt$tumor, header=T)
  colnames(tumor) = cutX(colnames(tumor))
  if (!startsWith(tumor[,1], "chr"))
    tumor[,1] <- paste0("chr", tumor[,1])
  tumor <- tumor[order(tumor[,1], as.numeric(tumor[,2])),]
  if (length(whichBedIsNA) > 0)
    tumor = tumor[-whichBedIsNA,]
  if (!checkBedAndCoverageValidity(bedFile, tumor)) {
    quit()
  }
  tumor <- as.matrix(tumor[,opt$colNum:ncol(tumor)])
  tumor = checkForDuplicatesAndRemove(tumor, opt$tumorSample)
  avgDepthTumorOn = determineAverageDepth(normal, bedFile)
  
  
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
  normalOff <- normalOff[-whichBedOffIsNA,]
  if (!checkBedAndCoverageValidity(bedFileOfftarget, normalOff)) {
    quit()
  }
  normalOff <- as.matrix(normalOff[,opt$colNum:ncol(normalOff)])
  # remain only samples that are in Normal cohort 
  normalOff <- checkForDuplicatesAndRemove(normalOff, opt$normalSample)
  
  
  normalOff <- normalOff[,which(colnames(normalOff) %in% colnames(normal))]
  avgDepthNormalOff = determineAverageDepth(normalOff, bedFileOfftarget)
  if (nrow(normalOff) != nrow(bedFileOfftarget)) {
    print("WARNING: your file with offtarget normal coverages have different amount of rows with offtarget bed file. It is most probably a technical mistake. Check the input.")
  }
  
  if (frameworkOff == "offtarget") {
    tumorOff <- ReadFileFast(opt$tumorOfftarget, header=T) 
    colnames(tumorOff) = cutX(colnames(tumorOff))
    if (!startsWith(tumorOff[,1], "chr"))
      tumorOff[,1] <- paste0("chr", tumorOff[,1])
    tumorOff <- tumorOff[order(tumorOff[,1], as.numeric(tumorOff[,2])),]
    # remain only samples that are in Tumor cohort 
    tumorOff <- tumorOff[-whichBedOffIsNA,]
    if (!checkBedAndCoverageValidity(bedFileOfftarget, tumorOff)) {
      quit()
    }
    tumorOff <- as.matrix(tumorOff[,opt$colNum:ncol(tumorOff)])
    tumorOff <- tumorOff[,which(colnames(tumorOff) %in% colnames(tumor))]
    avgDepthTumorOff = determineAverageDepth(tumorOff, bedFileOfftarget)
    
  }
}


print(paste("Started basic quality filtering.",Sys.time()))

rowsToRemove <- cleanDatasetFromLowCoveredFiles(normal, bedFile)
if (length(rowsToRemove) > 0) {
  toBind = cbind(bedFile[rowsToRemove,], rep("LowRawCoverage", length(rowsToRemove)))
  colnames(toBind)[ncol(toBind)] = "Description"
  bedPositionsThatWillBeFiltered = rbind(bedPositionsThatWillBeFiltered, toBind)
  
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

removeCentromeres = T
if (opt$hg38) {
  lstOfChromBorders <- getCytobands("cytobandsHG38.txt", removeCentromeres)
} else {
  lstOfChromBorders <- getCytobands("cytobands.txt", removeCentromeres)
}
left_borders <- lstOfChromBorders[[1]]
right_borders <- lstOfChromBorders[[2]]
ends_of_chroms <- lstOfChromBorders[[3]]


# check if any targets in BED are out of cytobands
for (chrom in unique(bedFile[,1])) {
  if (ends_of_chroms[[chrom]] < max(bedFile[bedFile[,1] == chrom,3])) {
    print("Coordinates in BED file are outside of the cytobands! Please check if your cytobands file matches your reference genome version!")
    quit()
  }
}

startX = NA
if (opt$par != "NO" & (framework == "germline" | frameworkOff == "offtargetGermline")) {
  modifiedListOfChromosomesWithPAR = addParalogousRegions(left_borders, right_borders, ends_of_chroms)
  startX = modifiedListOfChromosomesWithPAR[[1]]
  left_borders = modifiedListOfChromosomesWithPAR[[2]]
  right_borders = modifiedListOfChromosomesWithPAR[[3]]
  ends_of_chroms  = modifiedListOfChromosomesWithPAR[[4]]
}


if (frameworkDataTypes == "covdepthBAF" & !opt$onlyTumor) {
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
  overdispersionsNormal = listOfValues[[3]]
  overdispersionsTumor = listOfValues[[4]]
  pvaluesShifts = listOfValues[[5]]
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
if (nrow(lst[[3]]) > 0) {
  toBind = cbind(lst[[3]], rep("GCnormFailed", nrow(lst[[3]])))
  colnames(toBind)[ncol(toBind)] = "Description"
bedPositionsThatWillBeFiltered = rbind(bedPositionsThatWillBeFiltered, toBind)
}
normal <- lst[[1]]
if (framework == "somatic") {
  if (frameworkDataTypes == "covdepthBAF" & !opt$onlyTumor) {
    if (lengthBasedNorm)
      tumor <- lengthBasedNormalization(tumor, bedFile, allowedChroms=allowedChromsBaf)
    lst <- gc_and_sample_size_normalise(bedFile, tumor, allowedChroms=allowedChromsBaf)
  } else {
    tumor <- lengthBasedNormalization(tumor, bedFile)
    lst <- gc_and_sample_size_normalise(bedFile, tumor)
  }
  tumor <- lst[[1]]
  bedFile <- lst[[2]]
  writeOutLevelOfNoiseVersusCoverage(avgDepthTumorOn, tumor, bedFile, paste0(opt$out, "/ontargetTumor.summary.xls"))
}
bedFile <- lst[[2]]
writeOutLevelOfNoiseVersusCoverage(avgDepthNormalOn, normal, bedFile, paste0(opt$out, "/ontargetNormal.summary.xls"))


### OFF TARGET GC NORMALIZATION
if (frameworkOff == "offtarget" | frameworkOff == "offtargetGermline") {
  lst <- gc_and_sample_size_normalise(bedFileOfftarget, normalOff)
  normalOff <- lst[[1]]
  if (frameworkOff == "offtarget") {
    if (frameworkDataTypes == "covdepthBAF" & !opt$onlyTumor) {
      lst <- gc_and_sample_size_normalise(bedFileOfftarget, tumorOff, allowedChroms=allowedChromsBaf)
    } else {
      lst <- gc_and_sample_size_normalise(bedFileOfftarget, tumorOff)
    }
    tumorOff <- lst[[1]]
    bedFileOfftarget <- lst[[2]]
    writeOutLevelOfNoiseVersusCoverage(avgDepthTumorOff, tumorOff, bedFileOfftarget, paste0(opt$out, "/offtargetTumor.summary.xls"))
  }
  bedFileOfftarget <- lst[[2]]
  writeOutLevelOfNoiseVersusCoverage(avgDepthNormalOff, normalOff, bedFileOfftarget, paste0(opt$out, "/offtargetNormal.summary.xls"))
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
  toBind = cbind(bedFile[regionsToFilerOutOn,], rep("SystematicallyLowCov", length(regionsToFilerOutOn)))
  colnames(toBind)[ncol(toBind)] = "Description"
  bedPositionsThatWillBeFiltered = rbind(bedPositionsThatWillBeFiltered, toBind)
  normal = normal[-regionsToFilerOutOn,] + 10**-20
  if (framework == "somatic") {
    tumor = tumor[-regionsToFilerOutOn,] + 10**-20
  }
  bedFile = bedFile[-regionsToFilerOutOn,]
}
outputQCFailed = F

#print("These regions will be filtered out:")
#print(bedPositionsThatWillBeFiltered)

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
sdsForQC = apply(sqrt(normal[which(!bedFile[,1] %in% c("chrX","chrY")),]), 2, mad)
samplesToFilterOut = which(sdsForQC < 0.005 | sdsForQC > 0.5)
if (length(samplesToFilterOut) > 0) {
  print(paste("Germline samples", colnames(normal)[samplesToFilterOut], "did not pass QC due to high level of noise"))
  normal = normal[,-samplesToFilterOut]
}

print(paste("We start to cluster your data (you will find a plot if clustering is possible in your output directory)", opt$out, Sys.time()))
if (is.null(opt$clusterProvided)) {
  if (opt$noPlot) {
    clustering = rep("0", ncol(normal))
    names(clustering) = colnames(normal)
    print(clustering)
    outliersByClusteringCohort = c()
  } else {
    clusteringList <- returnClustering2(as.numeric(opt$minimumNumOfElemsInCluster))
    clustering = clusteringList[[1]]
    outliersByClusteringCohort = clusteringList[[2]]
  }
} else {
  clusteringList <- read.table(opt$clusterProvided, stringsAsFactors = F)
  clustering = clusteringList[,2]
  names(clustering) = clusteringList[,1]
  outliersByClusteringCohort = c()
}

orderOfBedFile <- order(bedFile[,1], as.numeric(bedFile[,2]))
bedFile = bedFile[orderOfBedFile,]
normal = normal[orderOfBedFile,]



print(paste("Gender estimation started", Sys.time()))
genderOfSamplesCohort <- Determine.gender(sqrt(normal), bedFile)

if (opt$sex != "" & !is.null(opt$normalSample)) {
  if (opt$normalSample %in% names(genderOfSamplesCohort)) {
    if (opt$sex %in% c("M","F")) {
      genderOfSamplesCohort[which(names(genderOfSamplesCohort) == opt$normalSample)] = opt$sex
    } else {
      print("--sex flag should be specified as M or F, and your flag is specified as")
      print(opt$sex)
      print("I quit.")
      quit()
    }
  } else {
    print("Your normalSample is not in the gender table! Thus, sex is not changed.")
  }
} else {
  if (opt$sex != "" & is.null(opt$normalSample)) {
    print("When you specify --sex flag, you need also to specify a sample name. I quit.")
    quit()
  }
}

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
    #if (cluster == -1) {
    #  print(paste("Samples from trio mode that are presented in trios.txt but do not have a full family in file", opt$normal , "will be excluded."))
    #  print(colnames(normal)[which(clustering == -1)])
    #  next
    #}
    
    

    
    
    samplesToAnalyse = which(clustering == cluster)
    if (!is.null(opt$normalSample)) {
      if (!opt$normalSample %in% colnames(normal)[samplesToAnalyse]) {
        next
      }
    }
    
    coverage <- sqrt(normal[,samplesToAnalyse])
    genderOfSamples = genderOfSamplesCohort[samplesToAnalyse]
    outliersByClustering = outliersByClusteringCohort[samplesToAnalyse]
    if (frameworkOff == "offtargetGermline") {
      samplesToAnalyseOff = which(colnames(normalOff) %in% colnames(coverage))
      coverageOff <- sqrt(as.matrix(normalOff))[,samplesToAnalyseOff]
      genderOfSamplesOff <- sapply(1:ncol(coverageOff), function(i) {return(genderOfSamplesCohort[which(colnames(normal) == colnames(coverageOff)[i])])})
    }
    
    
    
    #clusterExport(cl, c('EstimateModeSimple', 'bedFile', 'genderOfSamples', "lehmanHodges", 'Qn'))
    
    
    polymorphicRegions = NULL
    ### HERE THE POLYMORPHIC REGIONS DETECTION GOES
    if (opt$polymorphicCalling == "YES") {
      source(paste0(opt$folderWithScript, "/germline/mCNVsDetection.R"),local=TRUE)
      if (nrow(copyNumberForReportingGlobal) > 0)
        polymorphicRegions = copyNumberForReportingGlobal
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
    coverage.normalised = sweep(coverage, 1, mediansAndSds[[1]][,1] + 10**-40, FUN="/")
    rm(coverage)
    gc()
    
    
    
    sdsOfProbes = trimValues(as.numeric(mediansAndSds[[1]][,2]), 0.01)
    sdsOfGermlineSamples = mediansAndSds[[2]]
    filterOutRegionsWithSmallMedians <- which(mediansAndSds[[1]][,1] > 0.3)
    sdsOfProbes = sdsOfProbes[filterOutRegionsWithSmallMedians]
    bedFileFiltered = bedFile[filterOutRegionsWithSmallMedians,]
    coverage.normalised = coverage.normalised[filterOutRegionsWithSmallMedians,]
    print(paste("Amount of regions after low covered regions filtering (median < 0.3) in cluster", cluster, ":", round(100 * nrow(coverage.normalised) / numberOfRowsBeforeAllTheFiltrationNormal, digits = 3) ) )
    
    if (frameworkOff == "offtargetGermline") {
      
      autosomesOff = which(!bedFileOfftarget[,1] %in% c("chrX", "chrY"))
      mediansAndSdsOff = calculateLocationAndScale(bedFileOfftarget, coverageOff, genderOfSamplesOff, autosomesOff)
      coverage.normalised.off = sweep(coverageOff, 1, mediansAndSds[[1]][,1] + 10**-40, FUN="/")
      rm(coverageOff)
      gc()
      sdsOfProbesOff = trimValues(as.numeric(mediansAndSdsOff[[1]][,2]), 0.01)
      sdsOfGermlineSamplesOff = mediansAndSdsOff[[2]]
      filterOutRegionsWithSmallMedians <- which(mediansAndSdsOff[[1]][,1] > 0.3)
      sdsOfProbesOff = sdsOfProbesOff[filterOutRegionsWithSmallMedians]
      bedFileFilteredOfftarget = bedFileOfftarget[filterOutRegionsWithSmallMedians,]
      coverage.normalised.off = coverage.normalised.off[filterOutRegionsWithSmallMedians,]
      
    }
    
    if (!is.null(opt$triosFile)) {
      source(paste0(opt$folderWithScript,"/trios/germlineTrioSolver.R"),local=TRUE)
    } else if (opt$onlyTumor) {
      source(paste0(opt$folderWithScript,"/somatic/tumorOnlySolver.R"),local=TRUE)
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
  vect_of_norm_likeliks <- fast_dt_list(as.numeric(opt$degreesOfFreedomStudent))
  vect_of_t_likeliks <- fast_dt_list(as.numeric(opt$degreesOfFreedomStudent))
  
  samplesToAnalyse = which(clustering == cluster)
  genderOfSamples = genderOfSamplesCohort[samplesToAnalyse]
  tmpNormal = normal[,which(clustering == cluster)]
  if (frameworkOff == "offtarget") {
    tmpNormalOff = normalOff[,which(colnames(normalOff) %in% colnames(tmpNormal)),drop=F]
    bedFileForClusterOff = bedFileOfftarget
  }
  if (!is.null(opt$normalSample) & !is.null(opt$tumorSample)) {
    if (!opt$normalSample %in% colnames(tmpNormal)) {
      next
    }
  }
  bedFileForCluster = bedFile
  
  source(paste0(opt$folderWithScript, "/somatic/somaticSolver.R"),local=TRUE)
  if (frameworkOff == "offtarget") {
    prepareDataAndCall(bedFileForCluster, tmpNormal, tumor, genderOfSamples, bedFileForClusterOff, tmpNormalOff, tumorOff)
  } else {
    prepareDataAndCall(bedFileForCluster, tmpNormal, tumor, genderOfSamples)
  }
}


