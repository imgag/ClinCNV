
#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="folder with .cov files that have to be merged OR list of files with paths", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-n", "--numOfColumn"), type="integer", default=4, 
              help="column where the coverage starts [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input) | is.null(opt$out) | is.null(opt$numOfColumn)){
  print_help(opt_parser)
  stop("All arguments has to be specified in mergeFilesFromFolder.R script. Quit.", call.=FALSE)
}

merge_covs_from_folder <- function(folder, outputFileName, columnWhereCoveragesStart) {
  if (!file_test("-f", folder)) {
    filenames <- list.files(folder, pattern="*.cov", full.names=TRUE)
  } else {
    filenames <- as.vector(read.table(folder, header=F)[,1])
    indicesOfCovFilenames <- c()
    for (i in 1:length(filenames)) {
      if (substr(filenames[i], nchar(filenames[i]) - 3, nchar(filenames[i])) == ".cov") {
        indicesOfCovFilenames <- c(indicesOfCovFilenames, i)
      }
    }
    filenames = filenames[indicesOfCovFilenames]
  }
  
  firstRow <-read.table(filenames[1], nrows = 1, header = FALSE, sep ='\t', stringsAsFactors = FALSE, comment.char = "%")
  mergeBy = firstRow[1:(columnWhereCoveragesStart - 1)]
  mergeBy[1] = gsub("#", "X.", mergeBy[1])
  mergeBy = as.vector(unlist(mergeBy))
  
  mergeFilesIntoOne <- function(indices) {
    if (length(indices) > 1) {
      leftHalf <- indices[1:max(round(length(indices) / 2), 1)]
      rightHalf <- indices[(max(round(length(indices) / 2), 1) + 1):length(indices)]
      coverages1 <- mergeFilesIntoOne(leftHalf)
      coverages2 <- mergeFilesIntoOne(rightHalf)
      
      matrixOfMergedCoverages <-  merge(coverages1, coverages2,  by=mergeBy)
    } else if (length(indices) == 1) {
      matrixOfMergedCoverages = read.table(filenames[indices[1]], header=T,  comment.char = "%")
    } else {
      print("WTF")
    }
    
    return(matrixOfMergedCoverages)
  }
  
  coverages <- mergeFilesIntoOne(1:length(filenames))
  write.table(coverages, quote=F, sep="\t", row.names=F, file=outputFileName)
}


merge_covs_from_folder(opt$input, opt$out, opt$numOfColumn)
