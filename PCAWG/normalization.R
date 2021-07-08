setwd("/users/so/gdemidov/CNV/")
args = commandArgs(trailingOnly=TRUE)
chromosomesListPosition = as.numeric(args[1])
print(paste("Position in list:", chromosomesListPosition))

library(robustbase)
library(MASS)
library("data.table")
library(foreach)
library(doParallel)
no_cores <- min(detectCores() - 1, 4)
no_cores = 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)



cytobands <- read.table("/users/so/gdemidov/CNV/wholeWGS/cytobands.txt",stringsAsFactors = F, header = F, sep="\t", row.names=NULL)
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
gc_file <- read.table("/users/so/gdemidov/CNV/WGS/gc_one_kb_windows.txt",stringsAsFactors = F, header = F, sep="\t", row.names=NULL)

estimate_mode <- function(x) {
  d <- density(x, bw="SJ")
  #d <- density(x)
  d$x[which.max(d$y)]
}



paths <- c("/users/so/gdemidov/CNV/1st_batch_1_990", "/users/so/gdemidov/CNV/1st_batch_991_1800", "/users/so/gdemidov/CNV/1st_batch_1801_end", "/users/so/gdemidov/CNV/2nd_batch", "/users/so/gdemidov/CNV/3rd_batch")
listOfChromsomomes = c(1:22, "X", "Y")


# COMBINATIONS OF CHROMOSOMES FOR COMPUTATIONAL EFFECTIVENESS
chromosomesCombinationForMemoryLimit = list(c(1), c(2), c(3), c(4, 22), c(5, 21), c(6, 19), c(7, 20), c("X", 18), c(8, 17), c(10, 16), c(11, 15), c(12, 14), c(9, 13))
listOfChromsomomes = chromosomesCombinationForMemoryLimit[[chromosomesListPosition]]



#paths <- c(paths[1])

df <- "&"
for (i in 1:length(listOfChromsomomes)) {
  chr_df <- "&"
  for (path in paths) {
    setwd(path)
    print("Reading file...")
    chrom_full_name <- paste("chr", listOfChromsomomes[i], sep="")
    file_name <- paste("matrix_of_depths.", chrom_full_name, ".txt", sep="")
    print(file_name)
    colnames <- strsplit(readLines(file_name, n=1), "\t")[[1]]
    colnames <- c("mappability", colnames)
    setnames(local_df <- as.data.frame(fread(file_name, skip=1, header=F, stringsAsFactors = F, sep="\t")), colnames)
    
    
    gcs <- gc_file[gc_file[,1] == listOfChromsomomes[i],]
    
    coverages <- as.matrix(local_df[,4:ncol(local_df)]) 
    rowsToKeep <- which(as.numeric(local_df[,1]) > 500)
    print(length(rowsToKeep)/nrow(coverages))
    local_df <- local_df[rowsToKeep,]
    
    gcs <- gcs[gcs[,2] %in% local_df$start,]
    head(gcs)
    left_arm_start <- left_borders[chrom_full_name]
    right_arm_start <- right_borders[chrom_full_name]
    chromosome_vector <- rep(listOfChromsomomes[i], nrow(local_df))
    arms_vector <- rep(0, nrow(local_df))
    for (j in 1:nrow(local_df)) {
      if (local_df$end[j] < left_arm_start) {
        arms_vector[j] <- 1
      }
      if (local_df$end[j] > right_arm_start) {
        arms_vector[j] <- 2
      }
    }
    local_df <- cbind(gcs[,3], chromosome_vector, arms_vector, local_df)
    
    
    colnames(local_df) <- c("GC", "chr", "arm", "mapping", colnames(local_df)[-c(1,2,3,4)])
    if (chr_df == "&") {
      chr_df = local_df
    } else {
      coords <- intersect(chr_df[,"start"], local_df[,"start"])
      chr_df = chr_df[which(chr_df$start %in% coords),]
      local_df = local_df[which(local_df$start %in% coords),]
      colnames_chr_df <- colnames(chr_df)
      colnames_local_df <- colnames(local_df)
      unique_columns <- which(!(colnames_local_df %in% colnames_chr_df))
      print(colnames_local_df[which((colnames_local_df %in% colnames_chr_df))])
      chr_df <- cbind(chr_df, local_df[,-c(1:6)])
    }
  }
  if (df == "&") {
    df = chr_df
  } else {
    df <- rbind(df, chr_df)
  }
}


coverages <- as.matrix(df[,7:ncol(df)])
coverages <- (coverages + 10^-10) / 1000
info <- as.matrix(df[,1:6])

### DIAGNOSTICS OF DUPLICATES
list_of_duplicates_to_remove = read.table("/users/so/gdemidov/CNV/duplicates.pcawg")


rm(df)
rm(local_df)
rm(chr_df)
gc()


gcs <- round(as.numeric(info[,1]), digits = 2)
uniques_gcs <- unique(gcs)

coverages <- sqrt(coverages)

# sample size normalisation
setwd("/users/so/gdemidov/CNV/pcawg_last/files/summaries")
tablesOfNormalizationFactors <- c()
for (i in 1:5) {
  tablesOfNormalizationFactors = c(tablesOfNormalizationFactors, as.vector(read.table(paste0("sampleSizeNormFactor",i))[,1]))
}

for (i in 1:ncol(coverages)) {
  coverages[,i] <- (coverages[,i] / (tablesOfNormalizationFactors[i]))
}

# gc contennt normalization
gcsAllowed = seq(from=0.0, to=1.0, by=0.01)
tablesOfGcFactors <- list()
for (i in 1:5) {
  gcFactors = read.table(paste0("gcFactors",i))
  gcsAllowed = intersect(gcsAllowed, row.names(gcFactors))
  tablesOfGcFactors[[i]] = gcFactors
}

for (i in 1:5) {
  tablesOfGcFactors[[i]] = tablesOfGcFactors[[i]][which(row.names(tablesOfGcFactors[[i]]) %in% gcsAllowed),]
  tablesOfGcFactors[[i]] = tablesOfGcFactors[[i]][order(row.names(tablesOfGcFactors[[i]])),]
}

gcTableForWholeDataset = tablesOfGcFactors[[1]][which(row.names(tablesOfGcFactors[[i]]) %in% gcsAllowed),]#as.matrix(do.call("cbind", tablesOfGcFactors))
for (i in 2:5) {
  tablesOfGcFactors[[i]] = tablesOfGcFactors[[i]][which(row.names(tablesOfGcFactors[[i]]) %in% gcsAllowed),]
  gcTableForWholeDataset = cbind(gcTableForWholeDataset, tablesOfGcFactors[[i]])
}
gcTableForWholeDataset = as.matrix(gcTableForWholeDataset)

gcs <- round(as.numeric(info[,1]), digits = 2)
rowsToRemain = which(gcs %in% as.numeric(gcsAllowed))
coverages = coverages[rowsToRemain,]
info = info[rowsToRemain,]
gcs <- round(as.numeric(info[,1]), digits = 2)
uniques_gcs <- as.numeric(row.names(gcTableForWholeDataset))

# SWEEP GCS
for (i in 1:nrow(coverages)) {
  current_gc <- which(uniques_gcs == gcs[i])
  coverages[i,] <- coverages[i,] / gcTableForWholeDataset[current_gc,]
}




coverages = coverages[,-as.vector(unlist(list_of_duplicates_to_remove))]


clustering = unlist(read.table("/users/so/gdemidov/CNV/clustering.pcawg"))
clustersToRemain = which(table(clustering) >= 50)


matr_of_distances <- as.matrix(read.table("/users/so/gdemidov/CNV/matr_of_dist.pcawg"))
matr_of_distances = matr_of_distances[which(clustering %in% clustersToRemain), which(clustering %in% clustersToRemain)]

hc <- hclust(as.dist(1 - matr_of_distances^2), method="complete")
num_of_clusters = 5
clustering <- cutree(hc, k=num_of_clusters)
table(clustering)
plot(hc, 
     main="Dissimilarity = mad(vector of differences)", xlab="")
rect.hclust(hc, k=num_of_clusters, border="red")









plotDataFile <- function(x, functionType) {
  png("data.file.png", type = "cairo", width = 2500, height = 1600, pointsize=32)
  if (identical(functionType, plot)) {
    functionType(x, pch=19, col=rgb(0,0,0,0.1))
  } else {
    functionType(x, breaks=100)
  }
  dev.off()
}


library(MASS)
FindRobustMeanAndStandardDeviation <- function(x, modeEstimated = NA) {
  density_of_x <-  density(x, bw="SJ", kernel="gaussian")
  if (is.na(modeEstimated)) {
    mu = density_of_x$x[which.max(density_of_x$y)]
  } else {
    mu = modeEstimated
  }
  
  closest_to_mu <- which.min(abs(density_of_x$x - mu))
  which_are_bigger <- which(density_of_x$y > density_of_x$y[closest_to_mu])
  density_of_x <- as.data.frame(cbind(density_of_x$x, density_of_x$y))
  colnames(density_of_x) <- c("x","y")
  density_of_x[which_are_bigger,] <- density_of_x$y[closest_to_mu]
  EF = max(density_of_x$y)
  
  lower_bound = min(density_of_x$x)
  lower_bound_differs = F
  bounded_on_lower_copy_nuber = which(density_of_x$x < sqrt(mu ** 2 - 1/4))
  if (length(bounded_on_lower_copy_nuber) > 0) {
    start_to_the_left <- max(bounded_on_lower_copy_nuber)
    
    previous_value = density_of_x$y[start_to_the_left]
    for (i in seq(from = start_to_the_left, to=1, by=-1)) {
      AB =  density_of_x$y[i]
      if ((AB > previous_value + 10**-10) | AB < 10**-10) {
        lower_bound = density_of_x$x[i]
        lower_bound_differs = T
        break
      } else {
        previous_value = AB
      }
    }
  } 
  
  upper_bound = max(density_of_x$x)
  upper_bound_differs = F
  bounded_on_higher_copy_nuber = which(density_of_x$x > sqrt(mu ** 2 + 1/4))
  if (length(bounded_on_higher_copy_nuber) > 0) {
    start_to_the_right <- min(bounded_on_higher_copy_nuber)
    previous_value = density_of_x$y[start_to_the_right]
    for (i in seq(from = start_to_the_right, to=length(density_of_x$x), by=1)) {
      AB =  density_of_x$y[i]
      if ((AB > previous_value + 10**-10 | AB < 10**-10)) {
        upper_bound = density_of_x$x[i]
        upper_bound_differs = T
        break
      } else {
        previous_value = AB
      }
    }
  }
  
  if (upper_bound_differs & lower_bound_differs) {
    dtnorm0 <- function(X, mean, sd, log = TRUE) {dtnorm(X, mean, sd, lower_bound, upper_bound,
                                                         log)}
  } else if (!upper_bound_differs & !lower_bound_differs) {
    return(matrix(c(median(x), Qn(x)), nrow=1))
  } else if (upper_bound_differs) {
    dtnorm0 <- function(X, mean, sd, log = FALSE) {dtnorm(X, mean, sd, lower = -10**10, upper=upper_bound,
                                                          log)}
  } else {
    dtnorm0 <- function(X, mean, sd, log = FALSE) {dtnorm(X, mean, sd, lower=lower_bound, upper=10**10,
                                                          log)}
  }
  QnX <- Qn(x)
  data = x[which(x >= lower_bound & x <= upper_bound)]
  if (is.na(modeEstimated)) {
    result <- tryCatch({fitdistr(data, dtnorm0, start=list(mean=mean(data), sd=sd(data))); return(nres$estimate)}
                       , error = function(e) {return(matrix(c(mu, QnX), nrow=1))})
  } else {
    result <- tryCatch({fitdistr(data, dtnorm0, fix.arg=list(mean=modeEstimated), start=list(mean=modeEstimated, sd=sd(data))); return(nres$estimate)}
                       , error = function(e) {return(matrix(c(mu, QnX), nrow=1))})
  }
  
  # Sometimes we miss one cluster and that's cause to increase of standard deviation
  result[2] = min(QnX, result[2], sd(x))
  return(result)
}

NormalizeVarStabDataWithPeer <- function(normalized.coverages.for.parameters.estimation.normalized.by.mode) {
  set.seed(100)
  library(peer)
  library(qtl)
  model=PEER()
  n_factors = min(20, round(ncol(normalized.coverages.for.parameters.estimation.normalized.by.mode) / 10))
  PEER_setNk(model, n_factors)
  normalized.coverages.for.parameters.estimation.normalized.by.mode <- t(normalized.coverages.for.parameters.estimation.normalized.by.mode)
  PEER_setPhenoMean(model, as.matrix(normalized.coverages.for.parameters.estimation.normalized.by.mode))
  n_iterations = 50
  PEER_setNmax_iterations(model, n_iterations)
  PEER_update(model)
  residuals = PEER_getResiduals(model)
  differencesOfResiduals <- rep(0, ncol(normalized.coverages.for.parameters.estimation.normalized.by.mode))
  for (i in 1:ncol(residuals)) {
    differencesOfResiduals[i] = median(normalized.coverages.for.parameters.estimation.normalized.by.mode[,i] - residuals[,i])
  }
  residuals = sweep(residuals, 2, differencesOfResiduals, FUN="+")
  return(residuals)
}











genders <- unlist(read.table("/users/so/gdemidov/CNV/genders.pcawg.txt"))
most.represented.gender <- names(which.max(table(genders)))
which.are.males <- which(genders == "M")
which.are.females <- which(genders == "F")
which.is.most.popular <- which(genders == most.represented.gender)
which.is.not.most.popular <- which(genders != most.represented.gender)



setwd("/nfs/users2/so/gdemidov/CNV/pcawg_last/files/")
write.table(coverages, file=paste("coverages", chromosomesListPosition, sep="_"), quote=F)
write.table(info, file=paste("info", chromosomesListPosition, sep="_"), quote=F)
for (cluster in clustering) {
  if (!file.exists(paste("normalized.data.list_elem", chromosomesListPosition ,"cluster", cluster, sep="_"))) {
    coverages = as.matrix(read.table(paste("coverages", chromosomesListPosition, sep="_")))
    normalized.coverage.corrected.gc <- coverages[,which(clustering == cluster)]
    rm(coverages)
    
    normalized.coverage.corrected.gc.and.peer = t(NormalizeVarStabDataWithPeer(normalized.coverage.corrected.gc))
    
    modesAndSdsOfProbs <- matrix(unlist(sapply((1):nrow(normalized.coverage.corrected.gc.and.peer), function(i){
      if (i %% 1000 == 0) print(paste("We were able to process", i, "lines for now,", 
                                      round(100 * i/nrow(normalized.coverage.corrected.gc.and.peer), 3), "percents"))
      homozygous_deletions = which(normalized.coverage.corrected.gc.and.peer[i,] <= 0.2)
      coverage.of.probe = normalized.coverage.corrected.gc.and.peer[i,]
      
      
      if (info[i,2] == "X") {
        homozygous_deletions = union(homozygous_deletions, which.is.not.most.popular)
      }
      if (info[i,2] == "Y") {
        homozygous_deletions = union(homozygous_deletions, which.are.females)
      }
      # 25 points are ± enough for accurate estimation
      if (length(coverage.of.probe) - length(homozygous_deletions) < 25){ #0.25 * length(coverage.of.probe)) {
        robDev = c(0,0)
      } else {
        if (length(homozygous_deletions) > 0) {
          robDev = FindRobustMeanAndStandardDeviation(coverage.of.probe[-homozygous_deletions])
        } else { 
          robDev = FindRobustMeanAndStandardDeviation(coverage.of.probe)
        }
      }
      robDev
    })), nrow=nrow(normalized.coverage.corrected.gc.and.peer), byrow=T)
    
    homozygouslyDeletedReions = which(modesAndSdsOfProbs[,1] < 0.25 & modesAndSdsOfProbs > 2.5)
    #write.table(homozygouslyDeletedReions, file=paste0("excludedCoords", cluster, ".pcawg"), quote=F)
    #normalized.coverage.corrected.gc.and.peer = normalized.coverage.corrected.gc.and.peer[-homozygouslyDeletedReions,]
    #info = info[-homozygouslyDeletedReions,]
    #rm(normalized.coverage.corrected.gc)
    
    # CAN BE UNCOMMENTED FOR FURTHER PROCESSING
    #problematicStandardDeviationThreshold = quantile( modesAndSdsOfProbs[,2], 0.995)
    #problematicProbes <- which(modesAndSdsOfProbs[,2] > problematicStandardDeviationThreshold)
    #modesAndSdsOfProbs[problematicProbes,2] = problematicStandardDeviationThreshold
    
    # library(modeest)
    # EstimateModeSimple <- function(x) {
    #   tsybakov(x, kernel="biweight", alpha=0.75, dmp=TRUE)
    # }
    # for (coord in problematicProbes) {
    #   #print(modesAndSdsOfProbs[coord, 1])
    #   homozygous_deletions = c(which(normalized.coverage.corrected.gc.and.peer[coord,] <= 0.05))
    #   coverage.of.probe = normalized.coverage.corrected.gc.and.peer[coord,]
    #   if (info[i,2] == "X") {
    #     homozygous_deletions = union(homozygous_deletions, which.is.not.most.popular)
    #   }
    #   if (info[i,2] == "Y") {
    #     homozygous_deletions = union(homozygous_deletions, which.are.females)
    #   }
    #   if (length(homozygous_deletions) > 0) {
    #     modesAndSdsOfProbs[coord, 1] = EstimateModeSimple(coverage.of.probe[-homozygous_deletions])
    #   } else { 
    #     modesAndSdsOfProbs[coord, 1] = EstimateModeSimple(coverage.of.probe)
    #   }
    #   modesAndSdsOfProbs[coord, 2] = problematicStandardDeviationThreshold
    #   #print(modesAndSdsOfProbs[coord, 1])
    # }
    # 
    
    
    write.table(modesAndSdsOfProbs, file=paste("modesAndSdsOfProbs.list_elem",chromosomesListPosition ,"cluster", cluster, sep="_"), quote=F, row.names=F, col.names=F)
    save.image(paste("normalized.data.list_elem",chromosomesListPosition,"cluster", cluster, sep="_"))
    
    next  
    
    
    
    modesToDivide = modesAndSdsOfProbs[,1] 
    modesToDivide[homozygouslyDeletedReions, 1] = 1
    normalized.coverages.for.parameters.estimation.normalized.by.mode <- sweep(normalized.coverage.corrected.gc.and.peer, 1, modesToDivide[,1], FUN="/")
    
    
    
    
    normalized.coverages.for.parameters.estimation = normalized.coverages.for.parameters.estimation.normalized.by.mode - 1
    if (most.represented.gender == "M") {
      normalized.coverages.for.parameters.estimation[which(info[,2] == "X"),] =  normalized.coverages.for.parameters.estimation[which(info[,2] == "X"),] + (1 - sqrt(1/2))
    }
    normalized.coverages.for.parameters.estimation[which(info[,2] == "Y"),] =  normalized.coverages.for.parameters.estimation[which(info[,2] == "Y"),] + (1 - sqrt(1/2))
    
    
    
    
    # Caculate sample level of noise
    autosomes = which(!info[,2] %in% c("X","Y"))
    sample.standard.deviations <- apply(normalized.coverages.for.parameters.estimation[setdiff(autosomes,homozygouslyDeletedReions),], 2, Qn)
    normalized.coverages.for.parameters.estimation <- sweep(normalized.coverages.for.parameters.estimation, 2, sample.standard.deviations, FUN="/")
    
    
    info(logger, "7th - correction of SDs with sample specific standard deviations")
    modesAndSdsOfProbsLast <- matrix(unlist(sapply((1):nrow(normalized.coverages.for.parameters.estimation), function(i){
      if (i %% 100 == 0) print(paste("We were able to process", i, "lines for now,", 
                                     round(100 * i/nrow(normalized.coverage.corrected.gc), 3), "percents"))
      homozygous_deletions = which(normalized.coverage.corrected.gc.and.peer[i,] <= 0.25)
      coverage.of.probe = normalized.coverages.for.parameters.estimation[i,]
      
      if (info[i,2] == "X") {
        homozygous_deletions = union(homozygous_deletions, which.is.not.most.popular)
      }
      if (info[i,2] == "Y") {
        homozygous_deletions = union(homozygous_deletions, which.are.females)
      }
      if (length(coverage.of.probe) - length(homozygous_deletions) < 0.25 * length(coverage.of.probe)) {
        robDev = c(0,0)
      } else {
        if (length(homozygous_deletions) > 0) {
          robDev = FindRobustMeanAndStandardDeviation(coverage.of.probe[-homozygous_deletions])
        } else { 
          robDev = FindRobustMeanAndStandardDeviation(coverage.of.probe)
        }
      }
      robDev
    })), nrow=nrow(normalized.coverages.for.parameters.estimation), byrow=T)
    
    write.table(modesAndSdsOfProbs, file=paste("modesAndSdsOfProbs.list_elem",chromosomesListPosition ,"cluster", cluster, sep="_"), quote=F, rownames=F, colnames=F)
    write.table(modesAndSdsOfProbsLast, file=paste("modesAndSdsOfProbsLast.list_elem",chromosomesListPosition ,"cluster", cluster, sep="_"), quote=F, rownames=F, colnames=F)
    write.table(sample.standard.deviations, file=paste("sample.standard.list_elem",chromosomesListPosition ,"cluster", cluster, sep="_"), quote=F, rownames=F, colnames=F)
    
    save.image(paste("normalized.data.list_elem",chromosomesListPosition ,"cluster", cluster, sep="_"))
    
    
    alpha = 0.05 / nrow(modesAndSdsOfProbsLast)
    ratiosOfSds <- modesAndSdsOfProbs[,2] / modesAndSdsOfProbsLast[,2]
    medianOfRatios <- median(ratiosOfSds)
    sdsOfSdsRatio <- Sn(ratiosOfSds)
    lowerBoundOfSds <- medianOfRatios - qnorm(0.99) * sdsOfSdsRatio
    
    returnSdsForSampleAndProbe <- function(i, j) { # i - coord of probe, j - coord of sample
      baseSd = modesAndSdsOfProbsLast[i,2]
      if (ratiosOfSds[i] < lowerBoundOfSds) {
        baseSd = modesAndSdsOfProbsLast[i,2] / lowerBoundOfSds
      }
      return(baseSd * sample.standard.deviations[j])
    }
  }
}







