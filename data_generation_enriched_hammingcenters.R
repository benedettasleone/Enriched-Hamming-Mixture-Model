ham_enr_gen <- function(K, S, p, n, m, sigma_values, center_values, folder_path = NULL, save = TRUE, seed = 10091995) {
  
  set.seed(seed)
  
  n_per_outer <- round(n * S / sum(S))
  
  Y <- matrix(0, nrow = n, ncol = p)
  true_outer <- rep(1:K, times = n_per_outer)
  true_inner <- vector(length = n)
  
  sample_Yi <- function(c_vec, sigma_vec, m) {
    sapply(1:p, function(j) {
      mode <- c_vec[j]
      sigma <- sigma_vec[j]
      unnormalized <- rep(1 / exp(1 / sigma), m)
      unnormalized[mode] <- 1
      prob <- unnormalized / sum(unnormalized)
      sample(1:m, 1, prob = prob)
    })
  }
  
  inner_counter <- 1
  row_idx <- 1
  
  for (k in 1:K) {
    for (s in 1:S[k]) {
      
      c_ms <- rep(center_values[inner_counter], p)
      sigma_ms <- rep(sigma_values[inner_counter], p)
      n_s <- round(n_per_outer[k] / S[k])
      print(k)
      print(s)
      print(c_ms)
      for (i in 1:n_s) {
        Y[row_idx, ] <- sample_Yi(c_ms, sigma_ms, m)
        true_inner[row_idx] <- inner_counter
        row_idx <- row_idx + 1
      }
      inner_counter <- inner_counter + 1
    }
  }
  
  
  if (save) {
    write.csv(Y, file.path(folder_path, paste0("simulated_data_seed", seed, ".csv")), row.names = FALSE)
    write.csv(data.frame(outer = true_outer, inner = true_inner),
              file.path(folder_path, paste0("true_clusters_seed", seed, ".csv")), row.names = FALSE)
  }
  
  list(data = Y, groundTruth = true_outer, true_inner = true_inner)
}
