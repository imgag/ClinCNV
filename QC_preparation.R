library(robustbase)
library(MASS)
library("data.table")
library(foreach)
library(doParallel)
library(parallel)
no_cores <- min(detectCores(), 4)
cl<-makeCluster(no_cores)
registerDoParallel(cl)

setwd("/users/so/gdemidov/CNV/scripts_and_final_data")
args = commandArgs(trailingOnly=TRUE)
first_chrom = args[1]
last_chrom = args[2]
cluster_no = args[3]
image_name <- paste("chr", first_chrom, last_chrom, cluster_no, sep="_")
print(image_name)
load(image_name)
rm(coverages)
gc()
############
# Cleaning probes that did not work
############

which_has_na <- which(colSums(is.na(residuals)) > 0)
print(which_has_na)
if (length(which_has_na) > 0) {
  print(length(which_has_na))
  residuals <- residuals[,-which_has_na]
  info <- info[-which_has_na,]
  modes_vector <- modes_vector[-which_has_na]
}

tenth_ninety_quantile = sapply(1:ncol(residuals),function(i) {return(quantile(residuals[,i], c(0.1, 0.9)))})

less_than_one <- which(tenth_ninety_quantile[2,] < 1 | tenth_ninety_quantile[1,]  > 1)
checked_values <- less_than_one
print(dim(residuals))

############
# Part with sd estimation for all samples
############
sample_deviations <- apply(residuals, 1, mad)
print("Sample Sds calculated")
normalized_samples <- residuals - 1
for (i in 1:ncol(normalized_samples)) {
  if (length(which(normalized_samples[,i] > 0)) < 5) {
    normalized_samples[,i] = normalized_samples[,i] + 1
    shift_length = mean(normalized_samples[ which(normalized_samples[,i] > 0.5 & normalized_samples[,i] < 0.9) ,i])
    if (is.na(shift_length)) {
      shift_length = median(normalized_samples[,i])
    }
    normalized_samples[,i] = normalized_samples[,i] - shift_length
  }
}

normalized_samples <- sweep(normalized_samples, 1, sample_deviations, FUN="/")

print("Samples were normalised with Sample Standard Deviations...")

std_estimate <- function(m) {
  left_x = min(normalized_samples[,m]) 
  right_x = max(normalized_samples[,m])
  
  if (sum(is.na(normalized_samples[,m])) > 0) {print(normalized_samples[,m]); print(m);return(0)}
  dens <- density(normalized_samples[,m], bw="SJ", adjust=0.9)
  
  density_estimated_right = F  
  density_estimated_left = F
  border_shift <- 1 #1 - sqrt(3/4)
  x = dens$x
  y = dens$y
  if (length(which(x < border_shift)) > 1) {
    density_estimated_right = T
    starting_point <- max(x[which(x <  border_shift)])
    starting_y <- dens$y[which(dens$x == starting_point)]
    for (x in x[which(x > + border_shift)]) {
      current_y = dens$y[which(dens$x == x)]
      if (starting_y > current_y) {
        starting_y = current_y
      } else {
        right_x = x
        break
      }
    }
  } else {
    print(m)
    print("WARNING HERE")
  }
  
  x = dens$x
  y = dens$y
  if (length(which(x > -border_shift)) > 1) {
    density_estimated_left = T
    starting_point <- min(x[which(x > - border_shift)])
    starting_y <- dens$y[which(dens$x == starting_point)]
    if (length(dens$x[which(dens$x < - border_shift)]) - 1 > 1) {
      for (i in seq(from=length(dens$x[which(dens$x < - border_shift)]) - 1, to=2, by=-1)) {
        x = dens$x[which(dens$x <  - border_shift)][i]
        current_y = dens$y[which(dens$x == x)]
        if (starting_y > current_y) {
          starting_y = current_y
        } else {
          left_x = x
          break
        }
      }
    }
  } else {
    print(m)
    print("WARNING HERE")
  }
  dist_from_one <- max(- left_x, right_x)
  values_to_return <- which(normalized_samples[,m] >  - dist_from_one & normalized_samples[,m] <  dist_from_one)
  if (length(values_to_return) > 0.1 * nrow(residuals)) return(values_to_return)
  else {
    print(length(values_to_return))
    print(values_to_return)
    print(length(which(normalized_samples[,m] > -3 * dist_from_one & normalized_samples[,m] <  3 * dist_from_one)))
    return(which(normalized_samples[,m] > -3 * dist_from_one & normalized_samples[,m] <  3 * dist_from_one))
  }
  
}


probe_stand_devs <- rep(0, ncol(residuals))
probe_stand_devs_before <- rep(0, ncol(residuals))
probe_stds = matrix(0, ncol=0, nrow=2)


print("Calculated values for Standard Deviation calculation")
length_of_fragment = 10000
probe_stds = matrix(0, nrow=2, ncol=0)
for (j in seq(from=1, to=ncol(residuals), by=length_of_fragment)) {
  print(j)
  probe_stds_tmp <- matrix(0, nrow=2, ncol=0)
  for(i in j:min(j + length_of_fragment - 1, ncol(residuals))){
    values_to_use <- std_estimate(i)
    probe_stand_dev <- Qn(normalized_samples[,i][ values_to_use ])
    probe_stand_dev_before <- Qn(residuals[values_to_use, i])
    column_after_before_std <- matrix(c(probe_stand_dev, probe_stand_dev_before), nrow=2, ncol=1)
    probe_stds_tmp <- cbind(probe_stds_tmp, column_after_before_std)
  }
  gc()
  probe_stds = cbind(probe_stds, probe_stds_tmp)
  print("Finished up to...")
  print(min(j + length_of_fragment, ncol(residuals)))
}
print(dim(probe_stds))
print(dim(residuals))

print("Finished probe stds calculation...")
rm(list_of_values_for_sd)

probe_stand_devs <- probe_stds[1,]
probe_stand_devs_before <- probe_stds[2,]
print("Almost everything is finished except saving")
rm(normalized_samples)
gc()
get_sds_pr <- function(i) {
  return(sample_deviations * probe_stand_devs[i])
}

get_sds_sam <- function(i) {
  return(sample_deviations[i] * probe_stand_devs)
}

get_sd_pr_sam <- function(i, j){
  return(sample_deviations[i] * probe_stand_devs[j])
}

get_multipler_sam <- function(i) {
  return(median(probe_stand_devs_before / get_sds_sam(i)))
}

multipliers <- sapply(1:nrow(residuals), get_multipler_sam)
print(multipliers)
print("Finished multipliers")

############
# QC Control
############
setwd("/users/so/gdemidov/CNV/scripts_and_final_data/qc_prepared")
number_of_outliers <- sapply(1:nrow(residuals), function(i) {
  sds_sam <- 2.6 * get_sds_sam(i);
  return(length(which(residuals[i,] > 1 + sds_sam )) + length(which(residuals[i,] < 1 - sds_sam )) )
})
number_of_outliers <- number_of_outliers / ncol(residuals)

qc_drop_file_name <- paste("qc", cluster_no, sep="_")
if (first_chrom == 1){
  print("Finding outliers...")
  QC_drop_list <- which(number_of_outliers > 0.015)
  dups <- apply(residuals, 1, function(x) {return(length(which(x > 2)))})
  number_of_dels <- sapply(1:nrow(residuals), function(i) {
    sds_sam <- 2.6 * get_sds_sam(i);
    return(length(which(residuals[i,] < 1 - sds_sam )) )
  })
  number_of_dups <- sapply(1:nrow(residuals), function(i) {
    sds_sam <- 2.6 * get_sds_sam(i);
    return(length(which(residuals[i,] > 1 + sds_sam )))
  })
  number_of_outliers <- (number_of_dups + number_of_dels) / ncol(residuals)
  QC_drop_list <- which(number_of_outliers > 0.05 | dups > 10 | number_of_dups/number_of_dels > 5 | number_of_dels/number_of_dups > 5)
  print("Number of samples that did not pass QC...")
  print(length(QC_drop_list))
  save(QC_drop_list, file=qc_drop_file_name)
} else {
  load(qc_drop_file_name)
}
qc_prepared_file_name = paste("after_qc", first_chrom, last_chrom, cluster_no, sep="_")
save.image(qc_prepared_file_name)
stopCluster(cl)
