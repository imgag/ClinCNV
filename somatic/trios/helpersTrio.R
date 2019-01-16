
plotFoundCNVs <- function(found_CNVs, toyLogFoldChange, toyBedFile, outputFolder, chrom, cn_states, toySizesOfPointsFromLocalSds, plottingOfPNGs) {
  vector_of_states = cn_states
  cnvsToOutput <- matrix(0, nrow=0, ncol=11)
  if (nrow(found_CNVs) > 0) {
    for (s in 1:nrow(found_CNVs)) {
      priority = 0
      if (( 
        (vector_of_states[found_CNVs[s,4], 1] != 2)
        | (chrom %in% c("chrX", "chrY", "X", "Y") & (vector_of_states[found_CNVs[s,4], 1] != vector_of_states[found_CNVs[s,4], 2] &
                                                     vector_of_states[found_CNVs[s,4], 1] != vector_of_states[found_CNVs[s,4], 3]))
      )) priority = 1
      if (( 
        (vector_of_states[found_CNVs[s,4], 1] != 2 &
          (vector_of_states[found_CNVs[s,4], 1] != vector_of_states[found_CNVs[s,4], 2] &
           vector_of_states[found_CNVs[s,4], 1] != vector_of_states[found_CNVs[s,4], 3]))
          ) 
        )
        priority = 2
      #CNV_name <- paste(chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], "CN:", vector_of_states[found_CNVs[s,4]], "-2ln(loglik):", found_CNVs[s,1], 
      #                  ", number of regions:", found_CNVs[s,3] - found_CNVs[s,2] + 1 ,
      #                  ", length:", toyBedFile[found_CNVs[s,3],3] - toyBedFile[found_CNVs[s,2],2])
      #CNV_name_to_write <- paste(colnames(toyLogFoldChange)[sam_no],  chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], "CN",vector_of_states[found_CNVs[s,4]], sep="_")
      
      
      vectorOfGeneNames = c()
      genesThatHasToBeSeparated = unique(toyBedFile[found_CNVs[s,2]:found_CNVs[s,3],5])
      for (i in 1:length(genesThatHasToBeSeparated)) {
        vectorOfGeneNames = c(vectorOfGeneNames, unlist(strsplit(genesThatHasToBeSeparated[i], split=",")))
      }
      vectorOfGeneNamesTrimmed = c()
      for (elem in vectorOfGeneNames) {
        vectorOfGeneNamesTrimmed = c(vectorOfGeneNamesTrimmed,trimws(elem) )
      }
      annotationGenes <- paste(unique(vectorOfGeneNamesTrimmed), collapse=",")
      CNVtoOut <- matrix(c(chrom, toyBedFile[found_CNVs[s,2],2], toyBedFile[found_CNVs[s,3],3], vector_of_states[found_CNVs[s,4], 1], vector_of_states[found_CNVs[s,4], 2], vector_of_states[found_CNVs[s,4], 3],
                           priority,
                           -1 * found_CNVs[s,1], 
                           found_CNVs[s,3] - found_CNVs[s,2] + 1,
                           format((toyBedFile[found_CNVs[s,3],3] - toyBedFile[found_CNVs[s,2],2]) / 1000, nsmall=3),
                           annotationGenes), nrow=1)
      cnvsToOutput = rbind(cnvsToOutput, CNVtoOut)
      
      length_of_repr <- 1000
      
      
      st <- found_CNVs[s,2]
      fn <- found_CNVs[s,3]
      
      pr = plottingOfPNGs
      if (pr) {
        plot_st <- max(1,st - length_of_repr)
        plot_fn <- min(length(toyLogFoldChange), fn + length_of_repr)
        png(filename=paste0(outputFolder, "/", paste0(CNV_name_to_write, ".png")), type = "cairo", width = 1536, height = 640)
        #bitmap(filename=paste0(outputFolder, "/", paste0(CNV_name_to_write, ".png")) ,width = 1024, height = 640, units = "px" )
        
        
        
        plot(toyLogFoldChange[plot_st:plot_fn], main=CNV_name, ylab="Copy Number", xlab=(paste("CNV within Chromosome Arm" )),
             ylim=c(0, 3), cex=toySizesOfPointsFromLocalSds[plot_st:plot_fn], yaxt='n')
        
        axis(side = 2, at = sqrt(cn_states/2), labels = cn_states)
        abline(v=c(st - plot_st, st - plot_st + fn - st + 1), col="red")
        #text(x = (st - 2 * plot_st + fn) / 2, y = (max(toyLogFoldChange[st:fn]) + 0.1), pos=ifelse(found_CNVs[s,4] > 2, 3, 1), 
        #     labels=c(paste("Genes affected:", paste(unique(bedFile[st:fn, 5]), collapse=", " ))), cex=2,
        #     family = "mono")
        #legend(x = (st - 2 * plot_st + fn) / 2, y = (max(toyLogFoldChange[st:fn]) + 0.1), 
        #       legend=c(paste("Genes affected:", paste(unique(bedFile[st:fn, 5]), collapse=", " ))), cex=2)
        
        
        abline(h=sqrt(cn_states/2),lty=2,col=colours,lwd=3)
        
        points((st - plot_st):(st - plot_st + fn - st) + 1, toyLogFoldChange[st:fn],col="black", pch=21,bg=colours[found_CNVs[s,4]], cex=toySizesOfPointsFromLocalSds[found_CNVs[s,2]:found_CNVs[s,3]])
        
        
        ### EACH POINTS WITH DISTANCE > 10 MB WILL BE SEPARATED BY VERTICAL LINE
        distanceBetweenPoints = 10 ** 6
        for (i in 2:nrow(toyBedFile)) {
          if (toyBedFile[i,2] - toyBedFile[i - 1,2] > distanceBetweenPoints) {
            abline(v = i - 0.5 - plot_st, lty=2, col="grey")
          }
        }
        
        minimumDist <- min(sqrt(cn_states/2) - max(toyLogFoldChange[max(st - 50, plot_st):min(fn + 50, plot_fn)]))
        closestLine <- which(sqrt(cn_states/2) - max(toyLogFoldChange[max(st - 50, plot_st):min(fn + 50, plot_fn)]) == minimumDist)
        text(x = (st - 2 * plot_st + fn) / 2, y = (sqrt(cn_states/2)[closestLine] + 0.1), pos=3, 
             labels=c(paste("Genes affected:", paste(unique(bedFile[st:fn, 5]), collapse=", " ))), cex=1.5,
             family = "mono")
        dev.off()
      }
      
    }
  }
  return(cnvsToOutput)
}