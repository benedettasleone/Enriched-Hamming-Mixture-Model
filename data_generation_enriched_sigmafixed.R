ham_enr_gen = function(K, S, p, nj, m, sigma_values, center_values, folder_path = NULL, save = TRUE, seed = 10091995){
  
  n_per_outer <- rep(0,K)  
  inner_counter <- 1
  for(k in 1:K){
    for(s in 1:S[k]){
      n_per_outer[k] = n_per_outer[k] + nj[inner_counter]
      inner_counter = inner_counter + 1
    }
  }
  
  set.seed(seed)
  n = sum(nj)
  
  
  Y <- matrix(0, nrow = n, ncol = p)
  true_outer <- rep(1:K, times = n_per_outer)
  true_inner <- vector(length = n)
  
  row_idx <- 1
  inner_counter <- 1
  
  
  sample_Yi <- function(c_vec, sigma_vec, m) {
    probs <- sapply(1:p, function(j) {
      mode <- c_vec[j]
      sigma <- sigma_vec[j]
      unnormalized <- rep(1 / exp(1 / sigma), m)
      unnormalized[mode] <- 1
      prob <- unnormalized / sum(unnormalized)
      sample(1:m, 1, prob = prob)
    })
    return(probs)
  }
  
  for (k in 1:K) {
     sigma_m <- rep(sigma_values[k], p)
     for (s in 1:S[k]) {
      c_ms <- rep(center_values[inner_counter],p)
      n_s <- nj[inner_counter]  
      
      for (i in 1:n_s) {
        Y[row_idx, ] <- sample_Yi(c_ms, sigma_m, m)
        true_inner[row_idx] <- inner_counter
        row_idx <- row_idx + 1
      }
      inner_counter = inner_counter + 1
    }
  }

Y <- Y[1:(row_idx - 1), ]
true_outer <- true_outer[1:(row_idx - 1)]
true_inner <- true_inner[1:(row_idx - 1)]

if (save){
sim_data_filename <- file.path(folder_path, paste0("simulated_data_seed", seed, ".csv"))
true_clusters_filename <- file.path(folder_path, paste0("true_clusters_seed", seed, ".csv"))

write.csv(Y, sim_data_filename, row.names = FALSE)
write.csv(data.frame(outer = true_outer, inner = true_inner), true_clusters_filename, row.names = FALSE)
}

results             = list()
results$data        = Y
results$groundTruth = true_outer
#results$trueCent    = trueCent
results$true_inner  = true_inner
return(results)
}