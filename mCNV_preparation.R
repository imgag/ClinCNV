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
first_chrom = as.numeric(args[1])
last_chrom = as.numeric(args[2])
len = last_chrom - first_chrom
cluster_no = 5

setwd("/users/so/gdemidov/CNV/scripts_and_final_data/qc_prepared")
qc_prepared_file_name <- paste("after_qc", first_chrom, last_chrom, cluster_no, sep="_")

load(qc_prepared_file_name)

qc_drop_file_name <- paste("qc", cluster_no, sep="_")
load(qc_drop_file_name)



residuals <- residuals[-QC_drop_list,]
multipliers <- multipliers[-QC_drop_list]
sample_deviations <- sample_deviations[-QC_drop_list]
############
# Functions for mixture modelling
############

fast_norm_list <- function() {
  values <- seq(from = 0.0, to = 10000.0, by=1)
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



clustering_possible <- function(vec, maximum, minimum, sd_to_start) {
  if (vec[3] < minimum | tail(vec, n=3) > maximum) {
    return(FALSE)
  }
  return(TRUE)
}

std_estimate_around_one <- function(dist_to_work) {
  left_x = min(dist_to_work) 
  right_x = max(dist_to_work)
  dens <- density(dist_to_work, bw="SJ")
  border_shift <- 0.015 
  x = dens$x
  y = dens$y
  starting_point <- max(x[which(x < 1 + border_shift)])
  starting_y <- dens$y[which(dens$x == starting_point)]
  for (x in x[which(x > 1 + border_shift)]) {
    current_y = dens$y[which(dens$x == x)]
    if (starting_y > current_y) {
      starting_y = current_y
    } else {
      right_x = x
      break
    }
  }
  
  x = dens$x
  y = dens$y
  starting_point <- min(x[which(x > 1 - border_shift)])
  starting_y <- dens$y[which(dens$x == starting_point)]
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
  dist_from_one <- max(1 - left_x, right_x - 1)
  
  return(Qn(dist_to_work[which(dist_to_work > 1 - dist_from_one & dist_to_work < 1 + dist_from_one)]))
}



########
# Mixutre modelling
########

loc1 <- c(2)
loc12 <- c(1, 2)
loc23 <- c(2, 3)
loc13 <- c(1,2,3)
loc24 <- c(2,3,4)
loc14 <- c(1,2,3,4)
loc26 <- c(2,3,4,5,6)
loc28 <- c(2,3,4,5,6,7,8)
loc210 <- c(2,3,4,5,6,7,8,9,10)
loc212 <- c(2,3,4,5,6,7,8,9,10,11,12)

loc46 <- c(4,5,6)
loc48 <- c(4,5,6,7,8)
loc410 <- c(4,5,6,7,8,9,10)
loc68 <- c(6,7,8)
loc610 <- c(6,7,8,9,10)
loc16 <- c(1,2,3,4, 5, 6)
loc18 <- c(1,2,3,4, 5, 6, 7, 8)
loc110 <- c(1,2,3,4, 5, 6, 7, 8,9,10)


locations <- list(loc1, loc12, loc23, 
                  loc13, loc24, loc14, 
                  loc26, loc28, loc210, 
                  loc212, loc46,
                  loc48, loc410, loc68, 
                  loc16, loc18, loc110, loc610)


locations_low_mode <- c(1,2,3,4,5,6,7,8,9,10,11,16,17,18)
locations_medium_mode <- c(1,2,3,4,5,6,7,8,9,10,11,16,17,18)
locations_high_mode <- c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)



super_small_likelik = -10^20
return_res <- function(vect_of_values, k) {
  res <- rep(-2 * super_small_likelik, length(locations))
  maximum <- max(vect_of_values)
  minimum <- min(vect_of_values)
  if (shapiro.test(vect_of_values)$p > 0.01 ) {
    res[1] = 2 * super_small_likelik
  } else {
    #moda <- modes_vector_arm_chr[k]
    vect_of_allowed_locations <- seq(1:length(locations))
    max_divisor <- 8
    if (modes_vector_arm_chr[k] > sqrt(3/2)) {
      possible_divisors <- 5:8
    } else {
      possible_divisors <- 5
    }
    possible_divisors <- c(2, 3, 4, possible_divisors)
    if (maximum > sqrt(2)) {
      possible_divisors <- c(1, possible_divisors)
    }
    
    eps = 0.01
    
    values_to_check <- which(vect_of_values > 0.5)
    hom_del_values <- vect_of_values[!values_to_check]
    if (length(hom_del_values) > 10) {
      m <- mean(hom_del_values)
      if (length(hom_del_values) < 20) {
        s <- min(sd_to_start)
      } else {
        s <- Sn(hom_del_values)
      }
      values_to_check <- which(vect_of_values > max(0.45, m + 4 * s))
    } 
    sd_to_start <- std_estimate_around_one(vect_of_values[values_to_check]) * multipliers
    max_sd_to_start <- max(sd_to_start)
    
    for (j in (1:length(vect_of_allowed_locations))) {
      if(min(locations[[vect_of_allowed_locations[j]]]) <= max(possible_divisors)){
        possible_divisors_loc <- possible_divisors[which(possible_divisors >= min(locations[[vect_of_allowed_locations[j]]]))]
        if (length(possible_divisors_loc) == 0) {
          possible_divisors_loc = c(2)
        }
        location <- locations[[vect_of_allowed_locations[j]]]
        for (l in 1:length(possible_divisors_loc)) {
          
          divisor <- possible_divisors_loc[l]
          if (sqrt(min(location) / divisor) >= 0.45) {
            worth_to_check = T
            if (length(location) >= 6) {
              worth_to_check = clustering_possible(sqrt(location/divisor), maximum, minimum, sd_to_start)
            }
            if (worth_to_check) {
              tmp_list <- likelihood_t_mixture(sqrt(location/divisor), vect_of_values, values_to_check, sd_to_start, snp_coord)
              current_res <- -2 * tmp_list[[1]] + (length(location) + 1) * log(length(values_to_check))
              weight_of_current_clustering = tmp_list[[2]]
              n <- length(weight_of_current_clustering)
              
              
              if (res[vect_of_allowed_locations[j]] > current_res) {
                res[vect_of_allowed_locations[j]] <- current_res
              }
            } 
          }
        } 
      }
    }
  }
  return(res) 
}


rm(residuals_1)
rm(residuals_2)
gc()


matrices_of_likeliks <- list()
for (h in 0:(2 * len + 1)) {
  gc()
  chrom_to_study = as.numeric(first_chrom) + h %/% 2
  arm_to_study = h %% 2 + 1
  print("Chrom:")
  print(chrom_to_study)
  print("Arm:")
  print(arm_to_study)
  
  if (length(which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)) > 0) {
    resid <- residuals[,which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)]
    start_of_resid <- min(which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study))
    modes_vector_arm_chr <- modes_vector[which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study)]
    local_info <- info[which(info[,"chr"] == chrom_to_study & info[,"arm"] == arm_to_study),]
    
    print(Sys.time())
    matrix_of_likeliks <- foreach(i = 1:ncol(resid), 
                                  .combine = rbind, .packages='robustbase')  %dopar%  
      return_res(resid[,i], i)
    print(Sys.time())
    gc()
    matrices_of_likeliks[[10 * chrom_to_study + arm_to_study]] = matrix_of_likeliks
    matrices_of_likeliks_file_name <- paste("matrices", first_chrom, last_chrom, cluster_no, sep="_")
    setwd("/users/so/gdemidov/CNV/scripts_and_final_data/matrices_likeliks_prepared")
    save.image(matrices_of_likeliks_file_name)
  }
}

matrices_of_likeliks_file_name <- paste("matrices", first_chrom, last_chrom, cluster_no, sep="_")
setwd("/users/so/gdemidov/CNV/scripts_and_final_data/matrices_likeliks_prepared")
save.image(matrices_of_likeliks_file_name)
quit()