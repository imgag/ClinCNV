library(robustbase)
library(MASS)
library("data.table")


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




find_likeliks_for_different_models <- function(list_of_samples_names, chr) {
  path_to_snps = "/users/so/gdemidov/CNV/snps/likeliks/"
  lst_of_different_models_likeliks <- list()
  num_of_snps <- 10^10
  header <- unlist(strsplit(readLines(con="/users/so/gdemidov/CNV/snps/likeliks/snps_likeliks_2_1.txt", n=1), split="\t"))
  indices_to_keep <- which(header  %in% as.character(list_of_samples_names))
  indices_to_keep <- c(1, indices_to_keep)
  for (i in 1:6) {
    print(i)
    file_name <- paste("snps", "likeliks", chr, i, sep="_")
    file_name <- paste0(path_to_snps , file_name, ".txt")
    setnames(snps_table <- round(as.data.frame(fread(file_name, skip=1, header=F, stringsAsFactors = F, sep="\t"))), as.vector(unlist(header)))
    snps_table <- snps_table[,indices_to_keep]
    snps_table[is.na(snps_table)] = -10^20
    snps_table[,1] = snps_table[,1] * 1000 + 1
    lst_of_different_models_likeliks[[i+1]] = snps_table
    if (i == 1) {
      lst_of_different_models_likeliks[[1]] = snps_table
    }
  }
  return(lst_of_different_models_likeliks)
}



fast_norm_list <- function() {
  values <- seq(from = 0.0, to = 10000.0, by=1)
  vect_of_t_likeliks <- dnorm(values / 1000)
  return((vect_of_t_likeliks))
}
vect_of_t_likeliks <- fast_norm_list()
return_likelik <- function(x) {
  x = as.vector(x)
  x = round(abs(x * 1000)) + 1
  x = replace(x, which(x >= length(vect_of_t_likeliks)), length(vect_of_t_likeliks) - 1)
  return(vect_of_t_likeliks[x])
}



form_matrix_of_likeliks_one_sample <- function(i, j, k) {
  vector_of_values <- resid[k,]
  vector_of_states <- sqrt(c(0,1,2,3,4,5,6))
  vector_of_states <- vector_of_states / vector_of_states[3]
  matrix_of_BFs <- matrix(0, nrow=(j - i + 1), ncol=length(vector_of_states))
  start <- 1
  end <- j - i + 1
  if (i == 1 & j == ncol(resid)) {
    matrix_of_BFs = sapply(1:ncol(matrix_of_BFs), function(l) {
      sds <- get_sds_sam(k)[(start_of_resid + start + i - 2):(start_of_resid + end + i - 2)]
      value = return_likelik((vector_of_values - max(0.01, vector_of_states[l])) / sds ) / sds + 10^-100
      if (l == 1) {
        value = value + return_likelik((vector_of_values - 0.1) / sds ) / sds + 10^-100
        value = value + return_likelik((vector_of_values - 0.0001) / sds ) / sds + 10^-100
        value = value / 3
      }
      return(-2 * log(value))
    })
  } else {
    for (crd in start:end) {
      matrix_of_BFs[crd,] = sapply(1:ncol(matrix_of_BFs), function(l) {
        new_sd <- get_sd_pr_sam(k, start_of_resid + crd + i - 2)
        value = return_likelik((vector_of_values[crd + i - 1] - max(0.01, vector_of_states[l])) / new_sd) / new_sd
        if (l == 1) {
          value = value + return_likelik((vector_of_values[crd + i - 1] - 0.1) / sds ) / sds + 10^-100
          value = value + return_likelik((vector_of_values[crd + i - 1] - 0.0001) / sds ) / sds + 10^-100
          value = value / 3
        }
        return(-2 * log(value))
      })
    }
  }
  
  return(matrix_of_BFs)
}





find_final_state <- function(start, end, initial_state) {
  super_small_likelik = -10^20
  sweeped_matrix <- sweep(matrix_of_likeliks[(start + 1):(end - 1),,drop=F], 1, matrix_of_likeliks[(start + 1):(end - 1),initial_state])
  res_within <- apply(sweeped_matrix, 2, sum)
  cn_state_by_central_points <- which.min(res_within)
  
  only_central=F
  sweeped_matrix <- sweep(matrix_of_likeliks[(start):(end),,drop=F], 1, matrix_of_likeliks[(start):(end),initial_state])
  res <- apply(sweeped_matrix, 2, sum)
  likelihood_score_including_all_tiles <- min(res[cn_state_by_central_points], res_within[cn_state_by_central_points])
  if (likelihood_score_including_all_tiles == res_within[cn_state_by_central_points]) {
    only_central=T
    print("MIN HAPPENS")
    print(res[cn_state_by_central_points])
    print(res_within[cn_state_by_central_points])
    print(paste("Choosen", likelihood_score_including_all_tiles))
    print("")
  }
  
  if (!only_central) {
    sweeped_matrix <- sweep(matrix_of_likeliks_read_depth_only[(start):(end),,drop=F], 1, matrix_of_likeliks_read_depth_only[(start):(end),initial_state])
    res <- apply(sweeped_matrix, 2, sum)
    likelihood_score_all_tiles_read_depth_only <- res[cn_state_by_central_points]
    for_output <- apply((matrix_of_likeliks[(start + 1):(end - 1),,drop=F] / -2) * log10(e), 2, sum)
  } else {
    sweeped_matrix <- sweep(matrix_of_likeliks_read_depth_only[(start + 1):(end - 1),,drop=F], 1, matrix_of_likeliks_read_depth_only[(start + 1):(end - 1),initial_state])
    res <- apply(sweeped_matrix, 2, sum)
    likelihood_score_all_tiles_read_depth_only <- res[cn_state_by_central_points]
    for_output <- apply((matrix_of_likeliks[(start + 1):(end - 1),,drop=F] / -2) * log10(e), 2, sum)
  }
  if (which.max(for_output) != cn_state_by_central_points) {
    print(for_output)
    quit()
  }
  
  if (sum( is.infinite(for_output)) > 0) {
    print(start)
    print(end)
    print(sam_no)
    quit()
  }
  return(c(likelihood_score_including_all_tiles, cn_state_by_central_points, likelihood_score_all_tiles_read_depth_only, for_output))
}

maxSubArraySum <- function(x){
  bestSoFar = 0
  bestNow = 0
  bestStartIndexSoFar = -1
  bestStopIndexSoFar = -1
  bestStartIndexNow = -1
  for (i in 1:length(x)) {
    value = bestNow + x[i]
    if (value > 0) {
      if (bestNow == 0) {
        bestStartIndexNow = i
      }
      bestNow = value
    }
    else
      bestNow = 0
    
    if (bestNow > bestSoFar) {
      bestSoFar = bestNow
      bestStopIndexSoFar = i
      bestStartIndexSoFar = bestStartIndexNow
    }
  }
  return(c(bestSoFar, bestStartIndexSoFar, bestStopIndexSoFar))
}

find_one_CNV <- function(j, k, main_state, threshold, matrix_of_likeliks, min_CNV_len) {
  # sweeping the likelihoods
  subset_matrix <- -1 * matrix_of_likeliks[j:k,,drop=F]
  matrix_of_BFs <- sweep(subset_matrix, 1, subset_matrix[,main_state], FUN="-")
  
  coords_of_CNVs <- c(0,0,0,0)
  best_bf <- 0
  value <- threshold
  sequence_for_iteration = seq(1:ncol(matrix_of_BFs))
  sequence_for_iteration = sequence_for_iteration[-main_state]
  for (i in sequence_for_iteration) {
    res <- maxSubArraySum(matrix_of_BFs[,i])
    res[2] <- j + res[2] - 1
    res[3] <- j + res[3] - 1
    if (res[1] > max(value, best_bf)) {
      best_bf <- res[1]
      coords_of_CNVs <- c(res, i)
    }
  }
  if (coords_of_CNVs != c(0,0,0,0)) {
    return(coords_of_CNVs)
  }
  return(c(0,0,0,0))
}



find_all_CNVs <- function(minimum_length_of_CNV, threshold, price_per_tile, initial_state, matrix_of_likeliks) {
  vector_of_regions <- matrix(c(-10, 1, nrow(matrix_of_likeliks), initial_state), nrow=1, ncol=4)
  found_CNVs <- matrix(nrow=0, ncol=11)
  i = 1
  counter = 0
  while(i <= nrow(vector_of_regions)){
    current_region_to_look_for_CNVs = vector_of_regions[i,]
    start = current_region_to_look_for_CNVs[2]
    end = current_region_to_look_for_CNVs[3]
    allowed_length = max(3, minimum_length_of_CNV)
    flag_not_found_or_too_short = F
    if (end - start > allowed_length){
      found_CNV = find_one_CNV(start, end, current_region_to_look_for_CNVs[4], threshold, matrix_of_likeliks, minimum_length_of_CNV)
      if (found_CNV[4] == 0)
        flag_not_found_or_too_short = T
    } else {
      flag_not_found_or_too_short = T
    }
    if (flag_not_found_or_too_short) {
      # if we do not segment further we add CNV to the list - found CNV is not significant or 
      result_CNV <- current_region_to_look_for_CNVs
      if (current_region_to_look_for_CNVs[4] != initial_state & end - start >= minimum_length_of_CNV){
        bf_and_state <- find_final_state(start, end, 3)
        result_CNV[1] = bf_and_state[1]
        result_CNV[4] = bf_and_state[2]
        likelik_score_read_depth_only <- bf_and_state[3]
        result_CNV <- c(result_CNV, bf_and_state[4:length(bf_and_state)])
        if (bf_and_state[2] != initial_state & bf_and_state[1] < min(-threshold,  -price_per_tile * (end - start + 1) )) {
          if (likelik_score_read_depth_only < -threshold) {
            found_CNVs <- rbind(found_CNVs, result_CNV)
          }
        }
      }
    } else if (found_CNV[4] != 0 & found_CNV[3] - found_CNV[2] < minimum_length_of_CNV) {
      # found CNV is too short!
      evaluate_segment_further = T
      print(found_CNV)
      if (found_CNV[4] >= 2 & found_CNV[4] <= 4) {
        for (k in found_CNV[2]:found_CNV[3]) {
          matrix_of_likeliks[k,found_CNV[4]] = matrix_of_likeliks[k,current_region_to_look_for_CNVs[4]]
        }
      } else {
        if (found_CNV[4] == 1) {
          # replace borders as they are heterozygous deletions
          matrix_of_likeliks[max((found_CNV[2] - 1), start), 1] = matrix_of_likeliks[max((found_CNV[2] - 1), start), 2]
          matrix_of_likeliks[min((found_CNV[3] + 1), end), 1] = matrix_of_likeliks[min((found_CNV[3] + 1), end), 2]
        }
        if (found_CNV[4] > 4) {
          # replace borders as they are heterozygous duplications
          matrix_of_likeliks[max((found_CNV[2] - 1), start), found_CNV[4]] = matrix_of_likeliks[max((found_CNV[2] - 1), start), 4]
          matrix_of_likeliks[min((found_CNV[3] + 1), end), found_CNV[4]] = matrix_of_likeliks[min((found_CNV[3] + 1), end), 4]
        }
        found_CNV_with_changed = find_one_CNV(start, end, current_region_to_look_for_CNVs[4], threshold, matrix_of_likeliks, minimum_length_of_CNV)
        if (found_CNV_with_changed[3] == found_CNV[3] & found_CNV_with_changed[2] == found_CNV[2]) {
          # evaluate_segment_further = F
          for (k in max((found_CNV[2] - 1), start):min((found_CNV[3] + 1), end)) {
            matrix_of_likeliks[k,found_CNV[4]] = matrix_of_likeliks[k,current_region_to_look_for_CNVs[4]]
          }
        }
      }
      matrix_for_calculations <- -matrix_of_likeliks[start:end,] + matrix_of_likeliks[start:end, 3]
      if (!is.null(nrow(matrix_for_calculations))) {
        matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
      } 
      determined_state <- which.max(matrix_for_calculations)
      BF <- max(matrix_for_calculations)
      current_region_to_look_for_CNVs <- c(BF, start, end, determined_state)
      if (evaluate_segment_further)
        vector_of_regions <- rbind(vector_of_regions, current_region_to_look_for_CNVs)
    } else {
      # found CNV is significant and long
      start_of_CNV <- found_CNV[2]
      end_of_CNV <- found_CNV[3]
      new_state <- found_CNV[4]
      # segment found CNV itself
      if (end_of_CNV - start_of_CNV >= minimum_length_of_CNV) {
        vector_of_regions <- rbind(vector_of_regions, found_CNV)
      }
      if (start_of_CNV - start >= minimum_length_of_CNV) { # if left part is big enough to add! CAUTION!!!
        matrix_for_calculations <- -matrix_of_likeliks[start:(start_of_CNV - 1),] + matrix_of_likeliks[start:(start_of_CNV - 1), current_region_to_look_for_CNVs[4]]
        if (!is.null(nrow(matrix_for_calculations))) {
          matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
        } 
        determined_state <- which.max(matrix_for_calculations)
        BF <- max(matrix_for_calculations)
        left_part <- c(BF, start, start_of_CNV - 1, determined_state)
        vector_of_regions <- rbind(vector_of_regions, left_part)
      }
      if (end - end_of_CNV >= minimum_length_of_CNV) { # if right part is big enough to add! CAUTION!!!
        matrix_for_calculations <- -matrix_of_likeliks[start:(start_of_CNV - 1),] + matrix_of_likeliks[start:(start_of_CNV - 1), current_region_to_look_for_CNVs[4]]
        if (!is.null(nrow(matrix_for_calculations))) {
          matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
        } 
        determined_state <- which.max(matrix_for_calculations)
        BF <- max(matrix_for_calculations)
        right_part <- c(BF, end_of_CNV + 1, end, determined_state)
        vector_of_regions <- rbind(vector_of_regions, right_part)
      }
    }
    i = i + 1
  }
  return(found_CNVs)
}



colours = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414)]




sum_matrix_of_likeliks_depth_and_snps <- function(matrix_of_BFs, vect_of_coords_from_cov, sample_index ) {
  
  for (j in 1:ncol(matrix_of_BFs)) {
    matr_of_snp_likeliks = c(0, list_of_likeliks[[j]][,sample_index] / 10)
    matr_to_sum <- matrix(0, nrow=nrow(matrix_of_BFs), ncol=1)
    matr_to_sum = matr_to_sum + matr_of_snp_likeliks[vect_of_coords_from_cov[,2] + 1]
    matrix_of_BFs[,j] = matrix_of_BFs[,j] - matr_to_sum
  }
  return(matrix_of_BFs)
}
















len <- as.numeric(last_chrom) - as.numeric(first_chrom)
folder_name <- ("/users/so/gdemidov/CNV/scripts_and_final_data/rCNV_results_snp_9oct2017/")
super_small_likelik = -10^20
print("Started working on...")
print(date())

# SNPs part
sample_names = rownames(residuals)
dict_of_samples = read.table("/users/so/gdemidov/CNV/snps/merged_mappings.txt")
samples_presented_in_cluster = dict_of_samples[which(dict_of_samples[,1] %in% sample_names),2]



previous_chrom=-1
for (h in 0:(2 * len + 1)) {
  gc()
  chrom_to_study = first_chrom + h %/% 2
  arm_to_study = h %% 2 + 1
  
  #if (chrom_to_study > 11 | (chrom_to_study == 11 & arm_to_study == 2)) {
  write("CHROM AND ARM TO STUDY", stderr())
  write(chrom_to_study, stderr())
  write(arm_to_study, stderr())
  
  if (previous_chrom != chrom_to_study) {
    list_of_likeliks <- find_likeliks_for_different_models(samples_presented_in_cluster, chrom_to_study)
    previous_chrom = chrom_to_study
  }
  
  threshold <- 60
  price_per_tile <- 10
  initial_state <- 3
  minimum_length_of_CNV = 2 # final CNVs will be +1 length
  
  
  if (length(which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)) > 0){
    rst_fld <- paste(chrom_to_study, arm_to_study, sep="_")
    setwd(paste(folder_name, rst_fld, sep=""))
    resid <- residuals[,which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)]
    start_of_resid <- min(which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study))
    modes_vector_arm_chr <- modes_vector[which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)]
    local_info <- info[which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study),]
    
    mCNVs_folder <- ("/users/so/gdemidov/CNV/scripts_and_final_data/mCNV_results/")
    folder = paste(chrom_to_study, arm_to_study, sep="_")
    mCNVs_name <- paste("mCNVs", chrom_to_study, arm_to_study, sep="_")
    mCNVs_name <- paste(mCNVs_name, "txt", sep=".")
    mCNVs_to_exclude <- read.table(paste(mCNVs_folder, folder, "/", mCNVs_name, sep=""), stringsAsFactors=F, header=T)
    starts = sapply(1:nrow(mCNVs_to_exclude), function(i) {return(min(which(local_info[,5] >= mCNVs_to_exclude[i,3] - 10)) - 2)})
    ends = sapply(1:nrow(mCNVs_to_exclude), function(i) {return(max(which(local_info[,6] <= mCNVs_to_exclude[i,4] + 10)) + 2)})
    seq_to_exclude_before <- c()
    
    for (i in 1:length(starts)) {
      seq_to_exclude_before <- c(seq_to_exclude_before, starts[i]:ends[i])
    }
    seq_to_exclude <- c()
    for (elem in seq_to_exclude_before) {
      if (elem >= 1 & elem <= ncol(resid))
        seq_to_exclude <- c(seq_to_exclude, elem)
    }
    seq_to_exclude <- unique(seq_to_exclude)
    
    vec_of_mcnv_coords = c()
    
    total_no_of_CNVs <- 0
    
    
    vect_of_coords_from_cov <- matrix(0, ncol=2, nrow=ncol(resid))
    for (i in 1:ncol(resid)) {
      coord_in_resid <- as.numeric(local_info[i,5])
      vect_of_coords_from_cov[i,1] = coord_in_resid 
      coord_in_snps <- which(list_of_likeliks[[3]]$pos == coord_in_resid)
      if (length(coord_in_snps) == 1) {
        vect_of_coords_from_cov[i,2] = as.numeric(coord_in_snps) - 1
      } else {
        vect_of_coords_from_cov[i,2] = 0
        if (length(coord_in_snps) > 1) {
          print(coord_in_snps)
        }
      }
    }
    gc()
    if (as.numeric(cluster_no) > 0) {
      for (sam_no in 1:nrow(residuals)) {
        sample_name <- rownames(residuals)[sam_no]
        print(rownames(residuals)[sam_no])
        potential_names = -1
        if (length(which(dict_of_samples[,1] == sample_name)) > 0)
          potential_names = as.character(dict_of_samples[which(dict_of_samples[,1] == sample_name),2])
        print("Potential names")
        print(potential_names)
        sample_index = -1
        
        if (length(potential_names) > 1) {
          for (k in 1:length(potential_names)) {
            if (length(which(colnames(list_of_likeliks[[1]]) == potential_names[k])) == 1) {
              sample_index <- which(colnames(list_of_likeliks[[1]]) == potential_names[k])
            } else if (length(which(colnames(list_of_likeliks[[1]]) == potential_names[k])) > 1) {
              sample_index <- which(colnames(list_of_likeliks[[1]]) == potential_names[1])
            }
          }
        } else if (length(potential_names) == 1) {
          if (length(which(colnames(list_of_likeliks[[1]]) == potential_names[1])) > 0)
            sample_index <- which(colnames(list_of_likeliks[[1]]) == potential_names[1])
          if (length(sample_index) > 1) {
            sample_index <- sample_index[1]
          }
        }
        if (length(sample_index) > 1) sample_index = sample_index[1]
        
        dict_to_output = c()
        matrix_of_likeliks <- form_matrix_of_likeliks_one_sample(1, ncol(resid), sam_no)
        matrix_of_likeliks_read_depth_only <- matrix_of_likeliks
        if (sample_index > -1)
          matrix_of_likeliks <- sum_matrix_of_likeliks_depth_and_snps(matrix_of_likeliks, vect_of_coords_from_cov, sample_index )
        
        matrix_of_likeliks[seq_to_exclude,] = 0
        
        found_CNVs <- as.matrix(find_all_CNVs(minimum_length_of_CNV, threshold, price_per_tile, initial_state, matrix_of_likeliks))
        if (nrow(found_CNVs) > 1)
          found_CNVs <- found_CNVs[order(found_CNVs[,2]),]
        total_no_of_CNVs = total_no_of_CNVs + nrow(found_CNVs)
        if (nrow(found_CNVs) > 0) {
          for (s in 1:nrow(found_CNVs)) {
            CNV_name <- paste("chr", chrom_to_study, local_info[found_CNVs[s,2],5], local_info[found_CNVs[s,3],6], "CN:", found_CNVs[s,4] - 1, "-2ln(loglik):", found_CNVs[s,1])
            CNV_name_to_write <- paste(rownames(residuals)[sam_no], "chr", chrom_to_study, local_info[found_CNVs[s,2],5], local_info[found_CNVs[s,3],6], "CN", found_CNVs[s,4] - 1, sep="_")
            ten_based_likeliks <- paste(as.vector(found_CNVs[s,5:ncol(found_CNVs)]), collapse=";")
            print(ten_based_likeliks)
            CNV_string_to_write <- paste(rownames(residuals)[sam_no], chrom_to_study, local_info[found_CNVs[s,2],5], local_info[found_CNVs[s,3],6], found_CNVs[s,2], found_CNVs[s,3], found_CNVs[s,1], found_CNVs[s,4] - 1, ten_based_likeliks, sep="\t")
            dict_to_output <- c(dict_to_output, CNV_string_to_write)
            CNV_name_to_write <- paste(CNV_name_to_write, ".png", sep="")
            
            length_of_repr <- 1000
            which_mCNVs_in_window <- which(starts > max(1, found_CNVs[s,2] - length_of_repr) & starts < min(ncol(resid), found_CNVs[s,3] + length_of_repr) |
                                             which(ends > max(1, found_CNVs[s,2] - length_of_repr) & ends < min(ncol(resid), found_CNVs[s,3] + length_of_repr)))
            st <- found_CNVs[s,2]
            fn <- found_CNVs[s,3]
            
            pr = F
            if (pr) {
              png(filename=CNV_name_to_write, type = "cairo", width = 640, height = 640)
              
              plot_st <- max(1,st - length_of_repr)
              plot_fn <- min(ncol(resid), fn + length_of_repr)
              plot(resid[sam_no, plot_st:plot_fn]^2 * 2, main=CNV_name, ylab="Copy Number", xlab=(paste("Borders of CNV +/-",length_of_repr,"points" )),
                   ylim=c(0,max(resid[sam_no, plot_st:plot_fn]^2 * 2) + 0.02), cex=0.9)
              abline(v=c(st - plot_st, st - plot_st + fn - st), col="red")
              
              
              seq_of_cnvs <- c()
              if (length(which_mCNVs_in_window) > 0){
                for (i in 1:length(which_mCNVs_in_window)) {
                  for (l in starts[i]:ends[i]) {
                    seq_of_cnvs <- c(seq_of_cnvs, l)
                  }
                }
              }
              abline(v=c(seq_of_cnvs - plot_st, st - plot_st + seq_of_cnvs - st), col="gray52",lty=3)
              abline(h=0:8,lty=2,col=colours,lwd=3)
              points((st - plot_st):(st - plot_st + fn - st), resid[sam_no,st:fn]^2*2,col="black", pch=21,bg=colours[found_CNVs[s,4]])
              
              dev.off()
            }
            
            if (length(dict_to_output) > 0) {
              path = (paste(folder_name, rst_fld, sep=""))
              output_file_name <- paste(rownames(residuals)[sam_no], cluster_no, "txt", sep=".")
              output_file_name <- paste(path, output_file_name, sep="/")
              fileConn <- file(output_file_name)
              writeLines(dict_to_output, fileConn)
              close(fileConn)
            }
          }
        }
      }}
    
    print(paste("Total number:", nrow(residuals), "samples"))
    print(date())
  }
}




