library(robustbase)
library(data.table)

setwd("/users/so/gdemidov/CNV/scripts_and_final_data/qc_prepared")
args = commandArgs(trailingOnly=TRUE)
first_chrom = args[1]
last_chrom = args[2]
cluster_no = args[3]
qc_prepared_file_name <- paste("after_qc", first_chrom, last_chrom, cluster_no, sep="_")
write(qc_prepared_file_name, stderr())
print(qc_prepared_file_name)
print(date())

load(qc_prepared_file_name)

proxy_file_with_mCNVs_name = "/users/so/gdemidov/CNV/scripts_and_final_data/mCNV_results_genotyping/final_merge_5.txt"

proxy_file_with_mCNVs = read.table(proxy_file_with_mCNVs_name, stringsAsFactors=F, sep="\t", header=T)

qc_list = read.table("/users/so/gdemidov/CNV/scripts_and_final_data/mCNV_results_genotyping/qc_passed_samples.txt", stringsAsFactors=F)

residuals <- residuals[which(rownames(residuals) %in% as.vector(qc_list[,1])),]
multipliers <- multipliers[which(rownames(residuals) %in% as.vector(qc_list[,1]))]

##### MIXTURE MODELLING


std_estimate_around_one <- function(dist_to_work, border_shift, bw_adjust) {
  left_x = min(dist_to_work) 
  right_x = max(dist_to_work)
  dens <- density(dist_to_work, bw="SJ", adjust=bw_adjust) #1 - sqrt(3/4)
  x = dens$x
  y = dens$y
  starting_point <- max(x[which(x < 1 + border_shift)])
  starting_y <- dens$y[which(dens$x == starting_point)]
  if (length(which(x > 1 + border_shift)) > 2) {
    for (x in x[which(x > 1 + border_shift)]) {
      current_y = dens$y[which(dens$x == x)]
      if (starting_y > current_y) {
        starting_y = current_y
      } else {
        right_x = x
        break
      }
    }
  }
  
  x = dens$x
  y = dens$y
  starting_point <- min(x[which(x > 1 - border_shift)])
  starting_y <- dens$y[which(dens$x == starting_point)]
  if (length(which(x < 1 - border_shift)) > 2) {
    for (i in seq(from=length(dens$x[which(dens$x < 1 - border_shift)]) - 1, to=2, by=-1)) {
      x = dens$x[which(dens$x < 1 - border_shift)][i]
      current_y = dens$y[which(dens$x == x)]
      if (starting_y > current_y) {
        starting_y = current_y
      } else {
        left_x = x
        break
      }
    }
  }
  dist_from_one <- max(1 - left_x, right_x - 1)
  
  return(Qn(dist_to_work[which(dist_to_work > 1 - dist_from_one & dist_to_work < 1 + dist_from_one)]))
}


find_medians_of_values <- function(start, end) {
  medians_of_values <- apply(resid[,min(start+1, end):max(start, end-1), drop=F], 1, median)
  return(medians_of_values)
}

fast_norm_list <- function() {
  values <- seq(from = 0.0, to = 5000.0, by=1)
  vect_of_t_likeliks <- dnorm(values / 1000)
  return((vect_of_t_likeliks))
}

vect_of_t_likeliks <- fast_norm_list()
return_student_likelik <- function(x) {
  x = as.vector(x)
  x = round(abs(x * 1000)) + 1
  x = replace(x, which(x >= length(vect_of_t_likeliks)), length(vect_of_t_likeliks) - 1)
  return(vect_of_t_likeliks[x])
}

robustify_likeliks <- function(points_likeliks, percentage_robust) {
  vect_sum <- rowSums(points_likeliks)
  
  return(list(vect_sum))
}




find_final_state <- function(start, end, initial_state, found_likeliks_snps, weights_from_biggest_cohort, locations_to_test) {
  separation <- F
  medians_of_values <- find_medians_of_values(start, end)
  medians_of_values <- medians_of_values / estimate_mode(medians_of_values[which(medians_of_values > 0.5)])
  moda <- median(modes_vector_arm_chr[min(start+1, end):max(start, end-1)])
  
  if (length(which(medians_of_values > 15)) > 0.1 * nrow(residuals)) {
    moda = estimate_mode(medians_of_values[which(medians_of_values > 15)])
    if (moda > 0.5) {
      vect_of_values <- medians_of_values / moda
    }
  } else {
    vect_of_values <- medians_of_values
  }
  
  
  possible_divisors <- locations_to_test[which(locations_to_test > 0)]
  if (min(possible_divisors) > 8) {
    possible_divisors <- 1:3
  }
  
  possible_divisors <- possible_divisors[which(possible_divisors > 0)]
  border_shift = 0.05
  bw_adjust = 0.8
  if (min(locations_to_test) > 2) {
    border_shift=0.02
    bw_adjust=0.3
  }
  sd_to_start <- std_estimate_around_one(vect_of_values, border_shift, bw_adjust)
  
  values_to_check <- which(vect_of_values > 0.45)
  hom_del_values <- vect_of_values[-values_to_check]
  m = 0
  if (length(hom_del_values) > 10) {
    m <- mean(hom_del_values)
    if (length(hom_del_values) < 20) {
      s <- max(max(sd_to_start), 0.1)
    } else {
      s <- Sn(hom_del_values)
    }
    values_to_check <- which(vect_of_values > max(0.45, m + 4 * s))
    hom_del_values <- vect_of_values[-values_to_check]
  }
  
  separation <- F
  divisor=2
  location <- locations_to_test
  likeliks_for_divisors = rep(10^20, length(possible_divisors))
  for (i in 1:length(possible_divisors)) {
    copy_number_likeliks <- abs(sweep(rep.row(medians_of_values, length(location)), 1, sqrt(location / possible_divisors[i]), FUN="-"))
    copy_number_likeliks <- t(apply(copy_number_likeliks, 1, function(x){dnorm(x, sd = multipliers * sd_to_start)}))
    copy_number_likeliks[which(copy_number_likeliks < dnorm(-5))] <- dnorm(-5)
    copy_number_likeliks <- sweep(copy_number_likeliks, 1, weights_from_biggest_cohort, FUN="*")
    likeliks_for_divisors[i] = sum(log(apply(copy_number_likeliks, 2, sum)))
  }
  best_divisor <- possible_divisors[which.max(likeliks_for_divisors)]
  
  copy_number = rep(3, length(medians_of_values))
  
  location <- 0:14
  location[1] = 0.1
  best_div <- best_divisor
  if (length(hom_del_values) > 5) location[1] = (mean(hom_del_values))^2
  single_obs_weight <- sum(weights_from_biggest_cohort) / (length(medians_of_values)) 
  weights <- rep(1 * single_obs_weight, length(location))
  
  for (j in 1:length(locations_to_test)) {
    weights[locations_to_test[j] + 1] = weights[locations_to_test[j] + 1] + weights_from_biggest_cohort[j]
  }
  weights = weights / sum(weights)
  
  vector_of_sds <- rep(0, nrow(residuals))
  vector_of_sds = multipliers * sd_to_start
  
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  copy_number_likeliks <- abs(sweep(rep.row(vect_of_values, length(location)), 1, sqrt(location / best_div), FUN="-"))
  
  copy_number_likeliks <- t(apply(copy_number_likeliks, 1, function(x){dnorm(x, sd = vector_of_sds)}))
  
  
  copy_number_likeliks <- sweep(copy_number_likeliks, 1, weights, FUN="*")
  
  log_copy_number_likeliks <- log10(copy_number_likeliks + 10^-100)
  print(log_copy_number_likeliks)
  
  copy_number <- apply(copy_number_likeliks, 2, which.max)
  
  print(Sys.time())
  print("End copy number calc")
  
  return(list(min(likeliks_for_divisors), which.min(likeliks_for_divisors), separation, best_divisor, copy_number, log_copy_number_likeliks))
}


rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


len <- as.numeric(last_chrom) - as.numeric(first_chrom)
folder_name <- ("/users/so/gdemidov/CNV/scripts_and_final_data/mCNV_results_genotyping/")
super_small_likelik = -10^20
previous_chrom = -1 

for (h in 0:(2 * len + 1)) {
  chrom_to_study = first_chrom + h %/% 2
  if (previous_chrom != chrom_to_study) {
    previous_chrom = chrom_to_study
  }
  arm_to_study = h %% 2 + 1
  print(chrom_to_study)
  print(arm_to_study)
  if (length(which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)) > 0) {
    rst_fld <- paste(chrom_to_study, arm_to_study, sep="_")
    print(rst_fld)
    setwd(paste(folder_name, rst_fld, sep=""))
    resid <- residuals[,which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)]
    start_of_resid <- min(which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study))
    modes_vector_arm_chr <- modes_vector[which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)]
    local_info <- info[which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study),]
    
    max_coord <- max(local_info[,6])
    min_coord <- min(local_info[,5])
    
    cnvs_coords <- proxy_file_with_mCNVs[which(proxy_file_with_mCNVs[,1] == chrom_to_study & as.numeric(proxy_file_with_mCNVs[,2]) >=
                                                 min_coord & as.numeric(proxy_file_with_mCNVs[,3]) <= max_coord),]
    
    
    initial_state <- 1
    
    row_names <- ""
    for (i in 1:length(rownames(resid))) {
      row_names <- paste(row_names, rownames(resid)[i], sep="\t")
    }
    header <- paste("div", "chr", "name", "start", "end", "min_CN", "max_CN", "number_of_tiles", "start_tile", "end_tile", "best_divisor", row_names, sep="\t")
    dict_to_output = c(header)
    for (j in 1:nrow( cnvs_coords)) {
      print(j)
      start = as.numeric(cnvs_coords[j,2])
      end = as.numeric(cnvs_coords[j,3])
      
      
      rows_to_study <- which(as.numeric(local_info[,5]) >= start & as.numeric(local_info[,6]) <= end)
      if (length(rows_to_study) > 2) {
        copy_number = unlist(cnvs_coords[j,4:ncol(cnvs_coords)]) - 1
        cluster_weights = table(copy_number)
        locations_to_test <- as.numeric(names(cluster_weights))
        weights_from_biggest_cohort <- as.numeric(cluster_weights) / sum(as.numeric(cluster_weights))
        
        start = rows_to_study[1]
        end = rows_to_study[length(rows_to_study)]
        coords_of_good_mappable_regions <- round(local_info[rows_to_study, 5] / 1000)
        print(start)
        print(end)
        
        
        
        copy_number <- rep(2, nrow(resid))
        CNV_name_to_write <- paste("CNV", chrom_to_study, start, end, min(copy_number) - 1, max(copy_number) - 1, cluster_no, sep="_")
        string_to_output = ""
        
        best_divisor <- 2
        
        found_likeliks_snps = NULL
        calc_copy_number=T
        genotyped_cnv <- find_final_state(start, end, initial_state, found_likeliks_snps, weights_from_biggest_cohort, locations_to_test)
        best_divisor <- genotyped_cnv[[4]]
        copy_number <- genotyped_cnv[[5]]
        log_copy_number_likeliks <- genotyped_cnv[[6]]
        
        string_of_log_copy_number_likeliks <- apply(log_copy_number_likeliks, 2, paste, collapse=";")
        
        medians_of_values <- apply(resid[,min(start+1, end):max(start, end-1), drop=F], 1, median)
        moda <- median(modes_vector_arm_chr[min(start+1, end):max(start, end-1)])
        medians_of_values <- medians_of_values / estimate_mode(medians_of_values[which(medians_of_values > 0.5)])
        
        location <- 0:14
        location[1] = 0.0001
        colours = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414, 337)]
        
        
        trans_values = medians_of_values^2 * best_divisor
        long_enough = (end-start + 1 >= 3)
        if (long_enough) {
          CNV_name <- paste("chr", chrom_to_study, local_info[start,5], local_info[end,6], "CN, from:", min(copy_number) - 1, "to:", max(copy_number) - 1)
          CNV_name_to_write <- paste("CNV", chrom_to_study, local_info[start,5], local_info[end,6], min(copy_number) - 1, max(copy_number) - 1, cluster_no, sep="_")
          CN_string <- ""
          for (i in 1:length(copy_number)) {
            corresponding_string_to_sample <- paste(copy_number[i], string_of_log_copy_number_likeliks[i], sep=";")
            CN_string <- paste(CN_string, corresponding_string_to_sample, sep="\t")
          }
          string_to_output = paste("0", chrom_to_study, CNV_name_to_write, local_info[start,5], local_info[end,6], min(copy_number) - 1, max(copy_number) - 1, end-start+1, start, end, best_divisor, CN_string, sep="\t")
          dict_to_output = c(dict_to_output, string_to_output)
        } else {
          CNV_name <- paste("SHORT CNV", "chr", chrom_to_study, local_info[start,5], local_info[end,6], "CN, from:", min(copy_number) - 1, "to:", max(copy_number) - 1)
          CNV_name_to_write  <- paste("short_CNV", chrom_to_study, local_info[start,5], local_info[end,6], min(copy_number) - 1, max(copy_number) - 1, cluster_no, sep="_")
          CN_string <- ""
          for (i in 1:length(copy_number)) {
            corresponding_string_to_sample <- paste(copy_number[i], string_of_log_copy_number_likeliks[i], sep=";")
            CN_string <- paste(CN_string, corresponding_string_to_sample, sep="\t")
          }
          
        }
        
        
        
        png(filename=CNV_name_to_write, type = "cairo", width = 640, height = 640)
        
        plot(trans_values, ylim=c(0, max(trans_values) + 0.5), ylab="Copy Number", xlab="Samples from the largest cohort", main=CNV_name)
        abline(h=location)
        points(1:length(copy_number), trans_values, col="black", pch=21,bg=colours[copy_number])
        dev.off()
        
        
      } else {
        copy_number = rep(-1, nrow(resid))
        CNV_name <- paste("chr", chrom_to_study, as.numeric(cnvs_coords[j,2]), as.numeric(cnvs_coords[j,3]), "CN, from:", min(copy_number) - 1, "to:", max(copy_number) - 1)
        CNV_name_to_write <- paste("CNV", chrom_to_study, as.numeric(cnvs_coords[j,2]), as.numeric(cnvs_coords[j,3]), min(copy_number) - 1, max(copy_number) - 1, cluster_no, sep="_")
        CN_string <- ""
        for (i in 1:length(copy_number)) {
          CN_string <- paste(CN_string, copy_number[i], sep="\t")
        }
        string_to_output = paste("0", chrom_to_study, CNV_name_to_write, as.numeric(cnvs_coords[j,2]), as.numeric(cnvs_coords[j,3]), min(copy_number) - 1, max(copy_number) - 1, end-start+1, start, end, best_divisor, CN_string, sep="\t")
        dict_to_output = c(dict_to_output, string_to_output)
      }
    }
    CNVs_file_name <- paste("mCNVs", chrom_to_study, arm_to_study, cluster_no, sep="_")
    CNVs_file_name <- paste(CNVs_file_name, ".txt", sep="")
    fileConn<-file(CNVs_file_name)
    writeLines(dict_to_output, fileConn)
    close(fileConn)
  }
  
}

