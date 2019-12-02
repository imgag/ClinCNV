library(quantregForest)


#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--pathToClassifier"), type="character", default=NULL, 
              help="path to quantreg classifier, example = /path/to/randomForests.RObj", metavar="character"),
  make_option(c("--input"), type="character", 
              help="input fle path, example = /path/Sample_cnvs_clincnv.tsv", metavar="character"),
  make_option(c("--output"), type="character", 
              help="path to output file, example: /path/Sample_annotated_cnvs_clincnv.tsv", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

pathToClassifier = opt$pathToClassifier
print(pathToClassifier)
load(pathToClassifier)

inputFile = opt$input

con <- file(inputFile, "r")
forHeader = readLines(con)
close(con)
headerOnly = c()
for (line in forHeader) {
  if (substr(line, 1, 1) == "#") {
    headerOnly = c(headerOnly, line)
  }
}
headerOnly[length(headerOnly)] = paste0(headerOnly[length(headerOnly)], "\tFP_probability")
tableCNVs = read.table(inputFile, header=F, sep="\t")

estimateFDR <- function(tableOfCNVs) {
  FDR = rep(1, nrow(tableOfCNVs))
  AverageLoglikScore <- tableOfCNVs[,5]
  AverageLoglikPerTile = AverageLoglikScore / (0.001 * (tableOfCNVs[,3] - tableOfCNVs[,2] - 1))
  AmountWithinTheCohort = 1
  AverageQval = tableOfCNVs[,10]
  AverageNumOfMarkers = tableOfCNVs[,6]
  AverageAF = tableOfCNVs[,8]
  AverageLoglikPerTileCorrect = AverageLoglikScore / AverageNumOfMarkers
  minCN = tableOfCNVs[,4]
  maxCN = tableOfCNVs[,4]
  AmountOfCNVsPerSampleSingletons = rep(length(which(AverageLoglikScore > 5)), nrow(tableOfCNVs))
  data_for_rf = cbind(AverageLoglikScore,AverageLoglikPerTile,AmountWithinTheCohort,AverageQval,
                     AverageNumOfMarkers,AverageAF, AverageLoglikPerTileCorrect, minCN, maxCN,AmountOfCNVsPerSampleSingletons)
  colnames(data_for_rf) = c("AverageLoglikScore","AverageLoglikPerTile","AmountWithinTheCohort", "AverageQval","AverageNumOfMarkers","AverageAF", "AverageLoglikPerTileCorrect", "minCN", "maxCN","AmountOfCNVsPerSampleSingletons")

  fdrBorders = c(seq(from=0.5, to=0.999, by=0.005))
  conditionalQuantiles = matrix(0, nrow=0, ncol=length(fdrBorders))
  for (i in 1:nrow(tableOfCNVs)) {
    print(i)
    if (tableOfCNVs[i,4] > 2) {
      conditionalQuantilesOneCNV  <- predict(randomForests[["DUP"]],  data_for_rf[i,,drop=F], what=fdrBorders)
    } else {
      conditionalQuantilesOneCNV  <- predict(randomForests[["DEL"]],  data_for_rf[i,,drop=F], what=fdrBorders)
    }
    conditionalQuantiles = rbind(conditionalQuantiles, conditionalQuantilesOneCNV )
    
  }
  fdrsRes <- sapply(1:nrow(conditionalQuantiles), function(i) {res = 2 * (1 - fdrBorders[max(which(conditionalQuantiles[i,] <= 0.5))]); ifelse(is.na(res), 1, res)})
  
  return(fdrsRes)
}

annotatedOut = cbind(tableCNVs, estimateFDR(tableCNVs))
annotatedOutName = opt$output
con <- file(annotatedOutName, "w")
for (line in headerOnly) {
  writeLines(con=con, line)
}
close(con)
write.table(annotatedOut, file=annotatedOutName, quote = F, row.names=F, sep="\t", col.names = F, append = T)