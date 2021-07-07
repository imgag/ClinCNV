library(robustbase)
library(MASS)
library("data.table")
library("mixtools")

setwd("/users/so/gdemidov/CNV/scripts_and_final_data")
args = commandArgs(trailingOnly=TRUE)
first_chrom = args[1]
last_chrom = args[2]
cluster_no = 5







matrices_of_likeliks <- paste("matrices", first_chrom, last_chrom, cluster_no, sep="_")
setwd("/users/so/gdemidov/CNV/scripts_and_final_data/matrices_likeliks_prepared")
load(matrices_of_likeliks)

print("Everything was loaded")
########
# CNV calling
########



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

find_one_CNV <- function(j, k, main_state, threshold, matrix_of_likeliks) {
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

library(RColorBrewer)

likelihood_t_mixture <- function(vect, vect_of_values, values_to_check, sd_to_start, snp_coord) {
  data <- vect_of_values[values_to_check]
  percentage_robust = 0.05 * length(data) / length(vect_of_values)
  posterior_probs <- rep(1, length(vect_of_values))
  if (length(sd_to_start) > 1)
    sd_to_start <- sd_to_start[values_to_check]
  mu_mult <- 1
  if (length(vect) == 1){
    
    cluster_weights = rep(1/ length(vect), length(vect))
    
    points_likeliks <- (sapply(1:length(vect), 
                               function(i) {
                                 data_for_mean = (data - vect[i]) / sd_to_start; 
                                 return(cluster_weights[i] * 
                                          with(as.data.frame(data_for_mean), 
                                               return_student_likelik(data_for_mean)) / sd_to_start)
                               }))
    
    vect_sum <- rowSums(points_likeliks)
    lower_quantile <- min(quantile(vect_sum, percentage_robust))
    for_replacement <- min(vect_sum[which(vect_sum >= lower_quantile)])
    vect_sum <- replace(vect_sum, which(vect_sum <= lower_quantile), for_replacement)
    return(list(sum(log(vect_sum)), c(1), mean(sd_to_start)))
  }
  percentage_robust = 0 * length(data) / length(vect_of_values)
  delta_loglik <- 100
  eps = 1
  cluster_weights = rep(1/ length(vect), length(vect))
  points_likeliks <- (sapply(1:length(vect), 
                             function(i) {
                               data_for_mean = (data - vect[i]) / sd_to_start; 
                               return(cluster_weights[i] * 
                                        with(as.data.frame(data_for_mean), 
                                             return_student_likelik(data_for_mean)) / sd_to_start)
                             }))
  
  robust_res <- robustify_likeliks(points_likeliks, percentage_robust)
  vect_sum <- robust_res[[1]]
  outlier_values <- c()
  
  previous_loglik <- sum(log(vect_sum))
  weights = points_likeliks * (1 / vect_sum) / length(data)
  cluster_weights <- colSums(weights)
  
  data_sd <- data
  weigths_sd <- weights
  
  sd_counter = sqrt(sum(sapply(1:length(vect), 
                               function(i) {
                                 weigths_sd[,i] * (data_sd - vect[i])^2
                               })))
  sd_to_start = multipliers[values_to_check] * sd_counter
  counter = 0
  while(delta_loglik > eps) {
    points_likeliks <- (sapply(1:length(vect), 
                               function(i) {
                                 data_for_mean = (data - vect[i]) / sd_to_start; 
                                 return(cluster_weights[i] * 
                                          with(as.data.frame(data_for_mean), 
                                               return_student_likelik(data_for_mean)) / sd_to_start)
                               }))
    robust_res <- robustify_likeliks(points_likeliks, percentage_robust)
    vect_sum <- robust_res[[1]]
    outlier_values <- c()
    
    current_loglik <- sum(log(vect_sum))
    delta_loglik <- current_loglik - previous_loglik
    counter = counter + 1
    
    
    if (delta_loglik < eps) {
      return(list(previous_loglik, cluster_weights, sd_counter))
    }
    previous_loglik = current_loglik
    
    
    
    
    weights = points_likeliks * (1 / vect_sum) 
    
    vect_mu <- sapply(1:length(vect), function(i){return(sum(weights[,i] * data)/sum(weights[,i]))})
    vect_mu <- vect_mu + 10^(-10)
    
    weights <- weights / length(data)
    
    
    cluster_weights <- colSums(weights)
    bias <- sum(cluster_weights * (vect / vect_mu))
    data_sd <- data
    weigths_sd <- weights
    sd_counter = sqrt(sum(sapply(1:length(vect), 
                                 function(i) {
                                   weigths_sd[,i] * (data_sd - vect[i])^2
                                 })))
    
    sd_to_start = multipliers[values_to_check] * sd_counter
    if (!is.na(bias)) {
      final_mult <- mu_mult * bias
      if (final_mult < 1.01 & final_mult > 0.99 ) {
        mu_mult = final_mult
        vect <- vect * bias
      } else {
        if (final_mult > 1.01) {
          mu_mult = 1.01
          vect <- vect * bias
        }
        if (final_mult < 0.99) {
          mu_mult = 0.99
          vect <- vect * bias
        }
      }
    }
    
  }
  return(list(previous_loglik, cluster_weights, sd_counter, posterior_probs))
}

std_estimate_around_one <- function(dist_to_work) {
  left_x = min(dist_to_work) 
  right_x = max(dist_to_work)
  dens <- density(dist_to_work, bw="SJ")
  border_shift <- 0.05
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

find_final_state <- function(start, end, initial_state, calc_copy_number) {
  separation <- F
  res <- rep(-2 * super_small_likelik, length(locations))
  best_divisor <- rep(0, length(locations))
  best_sds <- rep(0, length(locations))
  best_weights <- list()
  super_small_likelik = -10^20
  medians_of_values <- find_medians_of_values(start, end)
  medians_of_values <- medians_of_values / estimate_mode(medians_of_values)
  moda <- median(modes_vector_arm_chr[min(start+1, end):max(start, end-1)])
  
  vect_of_values <- medians_of_values
  
  vect_of_allowed_locations <- 1:length(locations)
  
  
  possible_divisors_low_mode <- c(1,2,3)
  possible_divisors_medium_mode <- c(2,3,4,5)
  possible_divisors_high_mode <- c(2,3,4,5,6,7,8)
  possible_divisors <- possible_divisors_medium_mode
  if (moda < 1) {
    possible_divisors <- possible_divisors_low_mode
  }
  if (moda > sqrt(3/2)) {
    possible_divisors <- possible_divisors_high_mode
  }
  sd_to_start <- std_estimate_around_one(vect_of_values)
  
  values_to_check <- which(vect_of_values > 0.45)
  hom_del_values <- vect_of_values[!values_to_check]
  if (length(hom_del_values) > 10) {
    m <- mean(hom_del_values)
    if (length(hom_del_values) < 20) {
      s <- min(sd_to_start)
    } else {
      s <- Sn(hom_del_values)
    }
    values_to_check <- which(vect_of_values > max(0.45, m + 4 * s))
    hom_del_values <- vect_of_values[!values_to_check]
  }
  if (length(hom_del_values) > 10) {
    vect_of_allowed_locations = c(1,2,4,6,15,16,17)
  }
  
  separation <- F
  for (j in (1:length(vect_of_allowed_locations))) {
    location <- locations[[vect_of_allowed_locations[j] ]]
    if (min(location) < max(possible_divisors)){
      possible_divisors = min(location):max(possible_divisors)
      for (l in 1:length(possible_divisors)) {
        divisor <- possible_divisors[l]
        tmp_list <- 0
        if (sqrt(min(location, divisor) >= 0.5)) {
          tmp_list <- likelihood_t_mixture(sqrt(location/divisor), vect_of_values, values_to_check, sd_to_start, snp_coord)
          current_res <- -2 * tmp_list[[1]] + (length(location) + 2) * log(length(values_to_check))
          if (location[which.max(tmp_list[[2]])] != divisor) {
            current_res <- -2 * super_small_likelik
          }
          else {
            which_clust_signif <- which(tmp_list[[2]] > 10 / nrow(resid))
            first_sign <- which_clust_signif[1]
            last_sign <- tail(which_clust_signif,n=1)
            clust_weight <- min(tmp_list[[2]][first_sign], tmp_list[[2]][last_sign])
            which_clust_signif_inside <- which(tmp_list[[2]] > clust_weight - 1/nrow(resid))
            if (length(tmp_list[[2]]) > 2 & (length(which_clust_signif_inside) < last_sign - first_sign + 1 | 
                                             length(locations[[vect_of_allowed_locations[j] ]]) - length(which_clust_signif) > 2)) {
              current_res <- -2 * super_small_likelik
            } else {
              
              if (length(tmp_list[[2]]) > 1) {
                sds <- tmp_list[[3]]
                dists <- sapply(1:(length(location) - 1), function(i) {return(sqrt(location[i+1]) - sqrt(location[i]))})
                if (min(dists) > 4 * sds) {
                  separation <- T
                } else {
                  current_res <- -2 * super_small_likelik
                }
              }
            }
          }
          if (res[vect_of_allowed_locations[j]] > current_res) {
            res[vect_of_allowed_locations[j]] <- current_res
            best_divisor[vect_of_allowed_locations[j]] <- divisor
            best_weights[[vect_of_allowed_locations[j]]] = tmp_list[[2]]
            best_sds[vect_of_allowed_locations[j]] = tmp_list[[3]]
          }
        }
      }
    }
  }
  res = res - res[initial_state]
  
  copy_number = rep(2, length(medians_of_values))
  
  if (which.min(res) != 1 & calc_copy_number) {
    location <- 0:14
    location[1] = 0.0
    best_div <- best_divisor[which.min(res)]
    if (length(hom_del_values) > 5) location[1] = (mean(hom_del_values))^2 * best_div
    best_wt <- best_weights[[which.min(res)]] * ( length(values_to_check) / length(medians_of_values) ) 
    single_obs_weight <- sum(best_wt) / (length(values_to_check)) 
    best_sd <- best_sds[which.min(res)]
    states <- locations[[which.min(res)]]
    weights <- rep(0.1 * single_obs_weight, length(location))
    
    for (j in 1:length(states)) {
      weights[states[j] + 1] = weights[states[j] + 1] + best_wt[j]
    }
    weights[1] = weights[1] + single_obs_weight * (length(medians_of_values) - length(values_to_check))
    
    vector_of_sds <- rep(0, nrow(residuals))
    vector_of_sds = multipliers * best_sd
    
    rep.row<-function(x,n){
      matrix(rep(x,each=n),nrow=n)
    }
    copy_number_likeliks <- abs(sweep(rep.row(medians_of_values, length(location)), 1, sqrt(location / best_div), FUN="-"))
    
    copy_number_likeliks <- t(apply(copy_number_likeliks, 1, function(x){dnorm(x, sd = vector_of_sds)}))
    
    copy_number_likeliks <- sweep(copy_number_likeliks, 1, weights, FUN="*")
    
    copy_number <- apply(copy_number_likeliks, 2, which.max)
    
    
    print(Sys.time())
    print("End copy number calc")
    
  }
  return(list(min(res), which.min(res), separation, best_divisor[which.min(res)], copy_number))
}




find_all_CNVs <- function(minimum_length_of_CNV, threshold, initial_state, matrix_of_likeliks) {
  vector_of_regions <- matrix(c(-10, 1, nrow(matrix_of_likeliks), initial_state), nrow=1, ncol=4)
  found_CNVs <- matrix(nrow=0, ncol=5)
  potential_noisy_CNVs <- matrix(nrow=0, ncol=5)
  copy_numbers <- list()
  i = 1
  counter = 0
  found_CNVs_before_merging <- matrix(nrow=0, ncol=5)
  calc_copy_number = T
  while(i <= nrow(vector_of_regions)){
    current_region_to_look_for_CNVs = vector_of_regions[i,]
    start = current_region_to_look_for_CNVs[2]
    end = current_region_to_look_for_CNVs[3]
    found_CNV = find_one_CNV(start, end, current_region_to_look_for_CNVs[4], threshold, matrix_of_likeliks)
    if (found_CNV[4] == 0 | (end - start <= max(3, minimum_length_of_CNV))) {
      # if we do not segment further we add CNV to the list
      result_CNV <- current_region_to_look_for_CNVs
      if (current_region_to_look_for_CNVs[4] != initial_state ){
        bf_and_state <- find_final_state(start, end, initial_state, calc_copy_number)
        result_CNV[1] = bf_and_state[[1]]
        result_CNV[4] = bf_and_state[[2]]
        separation <- bf_and_state[[3]]
        divisor <- bf_and_state[[4]]
        copy_number <- bf_and_state[[5]]
        
        if (bf_and_state[2] != initial_state & bf_and_state[1] < -100) {
          found_CNVs_before_merging <- rbind(found_CNVs_before_merging, c(result_CNV, divisor) )
          copy_numbers[[start]] = copy_number
          #resid[,start:end] <- resid[,start:end] * sqrt(divisor/2)
        }
      }
    } else {
      start_of_CNV <- found_CNV[2]
      end_of_CNV <- found_CNV[3]
      new_state <- found_CNV[4]
      # segment found CNV itself
      vector_of_regions <- rbind(vector_of_regions, found_CNV)
      if (start_of_CNV - start >= 1) { # if left part is big enough to add! CAUTION!!!
        matrix_for_calculations <- -matrix_of_likeliks[start:(start_of_CNV - 1),] + matrix_of_likeliks[start:(start_of_CNV - 1), current_region_to_look_for_CNVs[4]]
        if (!is.null(nrow(matrix_for_calculations))) {
          matrix_for_calculations <- apply(matrix_for_calculations, 2, sum)
        } 
        determined_state <- which.max(matrix_for_calculations)
        BF <- max(matrix_for_calculations)
        left_part <- c(BF, start, start_of_CNV - 1, determined_state)
        vector_of_regions <- rbind(vector_of_regions, left_part)
      }
      if (end - end_of_CNV >= 1) { # if right part is big enough to add! CAUTION!!!
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
  print(vector_of_regions)
  
  calc_copy_number = T
  found_CNVs_before_merging <- found_CNVs_before_merging[order(found_CNVs_before_merging[,2]),]
  current_to_merge = c()
  
  found_CNVs_after_merging <- matrix(nrow=0, ncol=5)
  copy_numbers_final <- list()
  
  merged_cnv = c()
  for (i in 1:(nrow(found_CNVs_before_merging))) {
    medians_of_values_1 = find_medians_of_values(found_CNVs_before_merging[i,2], found_CNVs_before_merging[i,3])
    if (i < nrow(found_CNVs_before_merging)) {
      medians_of_values_2 = find_medians_of_values(found_CNVs_before_merging[i + 1,2], found_CNVs_before_merging[i + 1,3])
    } else {
      medians_of_values_2 = medians_of_values_1
    }
    is_close = sum(abs(copy_numbers[[found_CNVs_before_merging[i,2]]] - copy_numbers[[found_CNVs_before_merging[min(i + 1, nrow(found_CNVs_before_merging)),2]]]))
    is_close = is_close / nrow(residuals)
    correlatio = cor(medians_of_values_2, medians_of_values_1)
    distance = 0
    if (i < nrow(found_CNVs_before_merging)) distance = found_CNVs_before_merging[i + 1,2] - found_CNVs_before_merging[i,2]
    if (is_close < 0.05 | (correlatio > 0.9 & distance < 5)) {
      if (length(merged_cnv) == 0) {
        merged_cnv = c(found_CNVs_before_merging[i,2], found_CNVs_before_merging[min(i + 1, nrow(found_CNVs_before_merging)),3])
      } else {
        merged_cnv = c(merged_cnv[1], found_CNVs_before_merging[min(i + 1, nrow(found_CNVs_before_merging)),3])
      }
    } else {
      if (length(merged_cnv) == 0) {
        copy_numbers_final[[found_CNVs_before_merging[i,2]]] = copy_numbers[[found_CNVs_before_merging[i,2]]]
        found_CNVs_after_merging <- rbind(found_CNVs_after_merging, found_CNVs_before_merging[i,])
      } else {
        bf_and_state <- find_final_state(merged_cnv[1], merged_cnv[2], 1, T)
        result_CNV <- c(0, merged_cnv[1], merged_cnv[2] ,0)
        result_CNV[1] = bf_and_state[[1]]
        result_CNV[4] = bf_and_state[[2]]
        separation <- bf_and_state[[3]]
        divisor <- bf_and_state[[4]]
        
        found_CNVs_before_merging <- rbind(found_CNVs_before_merging, c(result_CNV, divisor) )
        copy_numbers_final[[start]] = bf_and_state[[5]]
        merged_cnv <- c()
      }
    }
    if (length(merged_cnv) != 0) {
      bf_and_state <- find_final_state(merged_cnv[1], merged_cnv[2], 1, T)
      result_CNV <- c(0, merged_cnv[1], merged_cnv[2] ,0)
      result_CNV[1] = bf_and_state[[1]]
      result_CNV[4] = bf_and_state[[2]]
      separation <- bf_and_state[[3]]
      divisor <- bf_and_state[[4]]
      
      found_CNVs_before_merging <- rbind(found_CNVs_before_merging, c(result_CNV, divisor) )
      copy_numbers_final[[start]] = bf_and_state[[5]]
      merged_cnv <- c()
    }
    
    
    
    
    
    print("")
  }
  
  
  return(list(found_CNVs_after_merging, potential_noisy_CNVs, copy_numbers_final))
}

len <- as.numeric(last_chrom) - as.numeric(first_chrom)
folder_name <- ("/users/so/gdemidov/CNV/scripts_and_final_data/mCNV_results_sd_corrected/")
for (h in 0:(2 * len + 1)) {
  chrom_to_study = first_chrom + h %/% 2
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
    matrix_of_likeliks = matrices_of_likeliks[[10 * chrom_to_study + arm_to_study]]
    
    minimum_length_of_CNV = 1
    threshold <- 500
    initial_state <- 1
    
    result <- find_all_CNVs(minimum_length_of_CNV, threshold, initial_state, matrix_of_likeliks)
    found_CNVs <- result[[1]]
    found_CNVs <- found_CNVs[order(found_CNVs[,2]),]
    potential_noisy_CNVs <- result[[2]]
    copy_numbers <- result[[3]]
    
    
    row_names <- ""
    for (i in 1:length(rownames(resid))) {
      row_names <- paste(row_names, rownames(resid)[i], sep="\t")
    }
    header <- paste("chr", "name", "start", "end", "min_CN", "max_CN", "number_of_tiles", "start_tile", "end_tile", "best_divisor", row_names, sep="\t")
    dict_to_output = c(header)
    for (j in 1:nrow(found_CNVs)) {
      
      start = found_CNVs[j,2]
      end = found_CNVs[j,3]
      
      copy_number <- copy_numbers[[start]]
      CNV_name_to_write <- paste("CNV", chrom_to_study, start, end, min(copy_number) - 1, max(copy_number) - 1, sep="_")
      string_to_output = ""
      
      best_divisor <- found_CNVs[j,5]
      medians_of_values <- apply(resid[,min(start+1, end):max(start, end-1), drop=F], 1, median)
      moda <- median(modes_vector_arm_chr[min(start+1, end):max(start, end-1)])
      medians_of_values <- medians_of_values / estimate_mode(medians_of_values)
      
      location <- 0:14
      location[1] = 0.0001
      colours = colors()[c(30, 114, 518, 148, 93, 456, 459, 552, 256, 652, 373, 68, 6, 600, 414, 337)]
      
      print(start)
      print(end)
      trans_values = medians_of_values^2 * best_divisor
      long_enough = (end-start + 1 >= 3)
      if (long_enough) {
        CNV_name <- paste("chr", chrom_to_study, local_info[found_CNVs[j,2],5], local_info[found_CNVs[j,3],6], "CN, from:", min(copy_number) - 1, "to:", max(copy_number) - 1)
        CNV_name_to_write <- paste("CNV", chrom_to_study, local_info[found_CNVs[j,2],5], local_info[found_CNVs[j,3],6], min(copy_number) - 1, max(copy_number) - 1, sep="_")
        CN_string <- ""
        for (i in 1:length(copy_number)) {
          CN_string <- paste(CN_string, copy_number[i], sep="\t")
        }
        string_to_output = paste(chrom_to_study, CNV_name_to_write, local_info[found_CNVs[j,2],5], local_info[found_CNVs[j,3],6], min(copy_number) - 1, max(copy_number) - 1, end-start+1, start, end, best_divisor, CN_string, sep="\t")
        dict_to_output = c(dict_to_output, string_to_output)
      } else {
        CNV_name <- paste("SHORT CNV", "chr", chrom_to_study, local_info[found_CNVs[j,2],5], local_info[found_CNVs[j,3],6], "CN, from:", min(copy_number) - 1, "to:", max(copy_number) - 1)
        CNV_name_to_write  <- paste("short_CNV", chrom_to_study, local_info[found_CNVs[j,2],5], local_info[found_CNVs[j,3],6], min(copy_number) - 1, max(copy_number) - 1, sep="_")
        CN_string <- ""
        for (i in 1:length(copy_number)) {
          CN_string <- paste(CN_string, copy_number[i] - 1, sep="\t")
        }
        string_to_output = paste(chrom_to_study, CNV_name_to_write, local_info[found_CNVs[j,2],5], local_info[found_CNVs[j,3],6], min(copy_number) - 1, max(copy_number) - 1, end-start+1, start, end, best_divisor, CN_string, sep="\t")
        dict_to_output = c(dict_to_output, string_to_output)
      }
      
      
      
      png(filename=CNV_name_to_write, type = "cairo", width = 640, height = 640)
      
      plot(trans_values, ylim=c(0, max(trans_values) + 0.5), ylab="Copy Number", xlab="Samples from the largest cohort", main=CNV_name)
      abline(h=location)
      points(1:length(copy_number), trans_values, col="black", pch=21,bg=colours[copy_number])
      dev.off()
      
      
      
    }
    CNVs_file_name <- paste("mCNVs", chrom_to_study, arm_to_study, sep="_")
    CNVs_file_name <- paste(CNVs_file_name, ".txt", sep="")
    potential_noisy_CNVs_name <- paste("noisy", chrom_to_study, arm_to_study, sep="_")
    potential_noisy_CNVs_name <- paste(potential_noisy_CNVs_name, ".txt", sep="")
    fileConn<-file(CNVs_file_name)
    writeLines(dict_to_output, fileConn)
    close(fileConn)
    write(t(potential_noisy_CNVs), file=potential_noisy_CNVs_name)
  }
  
}





