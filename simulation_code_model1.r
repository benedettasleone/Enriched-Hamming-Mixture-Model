#=========================================================================================
# SIMULATION STUDY 
#=========================================================================================

rm(list = ls())

#=========================================================================================
# LIBRARY IMPORTS
#=========================================================================================

# Core clustering libraries
library(mcclust)
library(mcclust.ext)

# Visualization libraries
library(ggplot2)
library(plotly)
library(htmlwidgets)

# Data manipulation and analysis
library(FactoMineR)
library(tidyverse)
library(readxl)
library(writexl)

# Custom functions and C++ code
source("RfunctionsAFP.R")
source("data_generation_enriched.R")
Rcpp::sourceCpp("enriched_hamming_centerfixed.cpp") 

#=========================================================================================
# SIMULATION PARAMETERS
#=========================================================================================

# Simulation settings
seed_set <- 1
burnin <- 5000
totiter <- 15000

# Model parameters
K <- 3          # outer clusters
Si <- c(3,3,3)  # inner clusters per outer
SS <- sum(Si)
mm <- 4         # modalities per variable
p <- 70         # number of variables
nj <- c(15,15,15,50,50,50,15,50,15) # number of observations

# Sigma values for each inner cluster
sigma_values <- rep(0.8, SS)


# Hyperparameters
gamma_tilde_m <- 0.2247071
gamma_tilde_s <- 0.2247071
gamma_mm <- 5
gamma_ss <- 3

#=========================================================================================
# INITIALIZATION
#=========================================================================================

cat('Starting simulation\n')

# Initialize data storage and results
data_list <- list()
results_summary <- data.frame(
  seed = numeric(),
  gamma_m = numeric(),
  gamma_s = numeric(),
  gamma_tilde_m = numeric(),
  gamma_tilde_s = numeric(),
  RI = numeric(),
  PosteriorModeClusters = numeric(),
  clusters_VI = numeric(),
  stringsAsFactors = FALSE
)

# Color map for plots
color_map <- c(
  "#FF9999", "#FF3333", "#990000",   
  "#99FF99", "#33CC33", "#006600",   
  "#99CCFF", "#3366FF", "#003399"    
)
names(color_map) <- as.character(1:9)

#log_factorial vector, used when sampling M and S
precomputed_log_factorial = rep(0,10100)
for (i in 2:10100) {
  precomputed_log_factorial[i] = precomputed_log_factorial[i-1] + log(i);
}

#=========================================================================================
# HELPER FUNCTIONS
#=========================================================================================

create_output_directory <- function(seed, p, sigma_values, totiter, burnin, nj, 
                                    gamma_m, gamma_s, gmm_m, gmm_s) {
  
  sigmas <- paste(sigma_values, collapse = "_")
  output_dir <- paste0("result_seed",seed,"_p",p,
                         "_sigma",sigmas,
                         "totiter",totiter+burnin,"_n",sum(nj),"_gamma_m",gamma_m,
                         "_gamma_s",gamma_s,"_lambda_m",gmm_m,"_lambda_s",gmm_s)
  folder_path <- paste0(getwd(),"/result_seed", seed,"_p",p,"_sigma",sigmas,"totiter",totiter+burnin,"_n",sum(nj),"_gamma_m",gamma_m,"_gamma_s",gamma_s,"_lambda_m",gmm_m,"_lambda_s",gmm_s)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  return(list(output_dir = output_dir, folder_path = folder_path))
}

create_mca_plot <- function(data, clusters, folder_path, seed, file_suffix = "", 
                            width = 8, height = 6, dpi = 300) {
  library(FactoMineR)
  library(ggplot2)
  
  dir.create(folder_path, showWarnings = FALSE, recursive = TRUE)
  
  Y_cat <- as.data.frame(data)
  Y_cat[] <- lapply(Y_cat, as.factor)
  
  mca_result <- MCA(Y_cat, ncp = 2, graph = FALSE)
  coords <- as.data.frame(mca_result$ind$coord)
  names(coords) <- make.names(names(coords))
  coords$cluster <- as.factor(clusters)
  
  p <- ggplot(coords, aes(x = Dim.1, y = Dim.2, color = cluster)) +
    geom_point(size = 2) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95")
    ) +
    labs(title = paste("MCA Plot - Seed", seed),
         x = "Dimension 1",
         y = "Dimension 2",
         color = "Cluster")
  
  output_file <- file.path(folder_path, paste0("mca_plot_seed", seed, file_suffix, ".png"))
  ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  
}

initialize_parameters <- function(data_list, seed, K, SS, p, sigma_values, nj, M_init = 8, S_init = 4) {
  
  n <- sum(nj)
  # Random initialization
  M <- M_init
  S <- rep(S_init, M)
  w_un <- rep(1, M)
  q_un <- rep(1, sum(S))
  s <- rep(0, n)
  m <- sample(M, n, replace = TRUE, prob = w_un / sum(w_un))
  
  # Initialize inner cluster assignments
  for (mm in unique(m)) {
    nm <- length(s[m == mm])
    if (mm == 1) {
      a <- 0
      b <- S[1]
    } else {
      a <- sum(S[1:(mm - 1)])
      b <- a + S[mm]
    }
    q_temp <- q_un[(a + 1):b]
    s[m == mm] <- sample(seq((a + 1), b), nm, replace = TRUE,
                         prob = q_temp / sum(q_temp))
  }
  
  sigma <- matrix(0.5, nrow = p, ncol = sum(S))
  center <- matrix(data = 1, nrow = p, ncol = M) 

  return(list(M = M, S = S, s = s, m = m, w_un = w_un, q_un = q_un, 
            sigma = sigma, center = center))
}

save_results <- function(post_burnin_chain, output_dir) {
  # Save main results
  write.table(post_burnin_chain$m, file.path(output_dir, "upper_clustering.csv"), 
              sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(post_burnin_chain$s, file.path(output_dir, "lower_clustering.csv"), 
              sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(post_burnin_chain$M, file.path(output_dir, "M.csv"), 
              sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(post_burnin_chain$k, file.path(output_dir, "K.csv"), 
              sep = ", ", row.names = FALSE, col.names = FALSE)
  
  # Save RDS files
  saveRDS(post_burnin_chain$center, file.path(output_dir, "centers.rds"))
  saveRDS(post_burnin_chain$zeta, file.path(output_dir, "zeta.rds"))
  saveRDS(post_burnin_chain$sigma, file.path(output_dir, "sigmas.rds"))
  saveRDS(post_burnin_chain$S, file.path(output_dir, "S.rds"))
  saveRDS(post_burnin_chain$S_m, file.path(output_dir, "S_m.rds"))
  saveRDS(post_burnin_chain$w_un, file.path(output_dir, "w_un.rds"))
  saveRDS(post_burnin_chain$q_un, file.path(output_dir, "q_un.rds"))
}

generate_diagnostic_plots <- function(output_dir, seed, post_burnin_chain) {
  
  # Number of clusters traceplot
  K_trace <- scan(file.path(output_dir, "K.csv"), sep = ",")
  df_K <- data.frame(Iteration = 1:length(K_trace), K = K_trace)
  
  p3 <- ggplot(df_K, aes(x = Iteration, y = K)) +
    geom_line(linewidth = 1) +
    theme_bw() +
    labs(title = paste("Traceplot of K seed", seed),
         x = "Iteration", y = "K")
  
  ggsave(filename = file.path(output_dir, "traceplot_number_of_clusters.png"), 
         plot = p3, width = 6, height = 4)
  
  # Posterior distribution of clusters
  post_total_cls <- table(K_trace) / length(K_trace)
  df <- data.frame(
    cluster_found = as.numeric(names(post_total_cls)),
    rel_freq = as.numeric(post_total_cls)
  )
  
  p1 <- ggplot(data = df, aes(x = factor(cluster_found), y = rel_freq)) +
    geom_col(fill = "black") +
    labs(x = "Clusters Found", y = "Relative Frequency",
         title = paste("Posterior distribution of clusters seed", seed)) +
    theme_bw() +
    scale_x_discrete(drop = FALSE)
  
  ggsave(filename = file.path(output_dir, "posterior_clusters_plot.png"), 
         plot = p1, width = 6, height = 4)
  
  # Sum of S traceplot
  S_list <- readRDS(file.path(output_dir, "S_m.rds"))
  S_sum <- sapply(S_list, sum)
  df_S <- data.frame(Iteration = 1:length(S_sum), SumS = S_sum)
  
  p2 <- ggplot(df_S, aes(x = Iteration, y = SumS)) +
    geom_line(linewidth = 1) +
    theme_bw() +
    labs(title = paste("Traceplot of sum(S) seed", seed),
         x = "Iteration", y = "sum(S)")
  
  ggsave(filename = file.path(output_dir, "trace_number_inner_components.png"), 
         plot = p2, width = 6, height = 4)
  
  # # Log-likelihood traceplot
  # df_loglik <- data.frame(Iteration = 1:length(post_burnin_chain$loglik), 
  #                        LogLik = post_burnin_chain$loglik)
  # 
  # p_loglik <- ggplot(df_loglik, aes(x = Iteration, y = LogLik)) +
  #   geom_line(linewidth = 1, color = "black") +
  #   theme_bw() +
  #   labs(title = paste("Traceplot of Log-Likelihood seed", seed),
  #        x = "Iteration", y = "Log-Likelihood")
  # 
  # ggsave(filename = file.path(output_dir, "traceplot_loglikelihood.png"), 
  #        plot = p_loglik, width = 6, height = 4)
  
  return(as.numeric(names(post_total_cls)[which.max(post_total_cls)]))
}

#=========================================================================================
# MAIN SIMULATION LOOP
#=========================================================================================

for (seed in seed_set) {
  for(gmm_m in gamma_tilde_m){
    for(gmm_s in gamma_tilde_s){
      for(gamma_m in gamma_mm){
        for(gamma_s in gamma_ss){
          
          cat('Simulation', seed, '\n')
          
          # Create output directory
          dirs <- create_output_directory(seed, p, sigma_values, totiter, burnin, nj, 
                                          gamma_m, gamma_s, gmm_m, gmm_s)
          output_dir <- dirs$output_dir
          folder_path <- dirs$folder_path
          
          # Generate data
          data_list[[seed]] <- ham_enr_gen(K, Si, p, nj, mm, sigma_values,
                                           folder_path, seed=seed)
          
          Y <- data_list[[seed]]$data
          clusters <- as.factor(data_list[[seed]]$true_inner)
          
          # MCA plot for true clusters
          create_mca_plot(Y, clusters, folder_path, seed)
          
          # Initialize parameters
          set.seed(seed)
          init_params <- initialize_parameters(data_list, seed, K, Si, p, sigma_values, nj)
          
          # Extract initialized parameters
          M <- init_params$M
          S <- init_params$S
          s <- init_params$s
          m <- init_params$m
          w_un <- init_params$w_un
          q_un <- init_params$q_un
          sigma <- init_params$sigma
          center <- init_params$center
          
          # Prepare for MCMC
          out <- rearrange_ms(m, s, S, M)
          m <- out$m
          s <- out$s
          stom <- out$stom
          S <- out$S
          k <- out$k
          h_m <- out$h_m
          
          # Set hyperparameters
          mu <- gamma_m  #hyper Poisson prior over M
          nu <- gamma_s  #hyper Poisson prior over S_m
          sig_m <- gmm_m #hyper gamma prior jumps out 
          sig_s <- gmm_s #hyper gamma prior jumps in
          
          # Attribute sizes and priors
          attrsize <- apply(data_list[[seed]]$data, 2, function(x) {length(table(x))})
          u <- v <- vector(length = p)
          
          u[attrsize == 3] <- 5.00
          v[attrsize == 3] <- 0.25
          u[attrsize == 4] <- 4.50
          v[attrsize == 4] <- 0.25
          u[attrsize == 5] <- 4.25
          v[attrsize == 5] <- 0.25
          
          Y_matrix <- as.matrix(Y)
          
          cat('Starting C++ MCMC... \n')
          start_time <- Sys.time()
          
          # Run C++ MCMC
          post_burnin_chain <- run_mcmc_cpp(
            Y = Y_matrix,
            m_init = m,
            s_init = s,
            w_un_init = w_un,
            q_un_init = q_un,
            center_init = center,
            sigma_init = sigma,
            M_init = M,
            S_init = S,
            stom = stom,
            attrsize = attrsize,
            u = u,
            v = v,
            burnin = burnin,
            totiter = totiter,
            mu = mu,
            nu = nu,
            sig_m = sig_m,
            sig_s = sig_s,
            precomputed_log_factorial = precomputed_log_factorial,
            verbose = TRUE
          )
          
          end_time <- Sys.time()
          cat("C++ MCMC completed in", round(end_time - start_time, 3), "seconds\n")
          
          # Save results
          save_results(post_burnin_chain, output_dir)
          
          # Diagnostic plots
          post_mode_clusters <- generate_diagnostic_plots(output_dir, seed, post_burnin_chain)
          
          # Posterior similarity matrix and clustering
          C <- as.matrix(read.table(file.path(output_dir, "upper_clustering.csv"), 
                                    sep = ",", header = FALSE))
          psm <- comp.psm(C)
          VI <- minVI(psm)
          
          cat("Cluster Sizes:\n")
          print(table(VI$cl))
          cat("\nAdjusted Rand Index:", arandi(VI$cl, data_list[[seed]]$groundTruth, adjust = FALSE), "\n")
          
          # Save PSM plot
          png(filename = file.path(output_dir, "psm_plot.png"), width = 800, height = 600)
          myplotpsm2(psm, classes = data_list[[seed]]$groundTruth, ax = FALSE, ay = FALSE)
          dev.off()
          
          # #cRAMER-RAO dependence structure
          # CV = round(CramerV(data_list[[seed]]$data),4)
          # png(filename = file.path(output_dir, "cramer_complete.png"), width = 800, height = 600)
          # myplotpsm2(CV)#,zlim=c(0,1))           
          # dev.off() 
          # 
          # #matrices for each true cluster
          # cluster_ids <- unique(data_list[[seed]]$groundTruth) 
          # 
          # for (i in seq_along(cluster_ids)) {
          #   #seq_along(x) è preferibile perché funziona anche se x è vuoto (length(x) = 0) e non dà errori come 1:0 se scrivessi 1:len(x)
          #   cluster_label <- cluster_ids[i]
          #   
          #   CV_c <- round(CramerV(data_list[[seed]]$data[data_list[[seed]]$groundTruth == cluster_label, ]), 4)
          #   
          #   filename <- file.path(output_dir, paste0("cramer_cluster_true_", i, ".png"))
          #   png(filename = filename, width = 800, height = 600)
          #   
          #   myplotpsm2(CV_c)#, zlim = c(0, 1), main = paste("Cramér Cluster", i))  # se la funzione supporta 'main'
          #   
          #   dev.off()
          # }
          # 
          # #matrices for each cluster found
          cluster_ids <- unique(VI$cl)
          # 
          # for (i in seq_along(cluster_ids)) {
          #   #seq_along(x) è preferibile perché funziona anche se x è vuoto (length(x) = 0) e non dà errori come 1:0 se scrivessi 1:len(x)
          #   cluster_label <- cluster_ids[i]
          #   
          #   CV_c <- round(CramerV(data_list[[seed]]$data[data_list[[seed]]$groundTruth == cluster_label, ]), 4)
          #   
          #   filename <- file.path(output_dir, paste0("cramer_cluster_estimated_", i, ".png"))
          #   png(filename = filename, width = 800, height = 600)
          #   
          #   myplotpsm2(CV_c)#, zlim = c(0, 1), main = paste("Cramér Cluster", i))  # se la funzione supporta 'main'
          #   
          #   dev.off()
          # }
          
          # Create clustered data plot
          create_mca_plot(Y, as.factor(VI$cl), folder_path, seed, "_clustered")
          
          # Save data and ground truth
          write_xlsx(as.data.frame(data_list[[seed]]$data), 
                     file.path(folder_path, paste0("dati_simulatiseed", seed, "_data.xlsx")))
          write_xlsx(as.data.frame(data_list[[seed]]$groundTruth), 
                     file.path(folder_path, paste0("dati_simulatiseed", seed, "_GT.xlsx")))
          
          # Update results summary
          ari_value <- arandi(VI$cl, data_list[[seed]]$groundTruth, adjust = FALSE)
          results_summary <- rbind(results_summary, data.frame(
            seed = seed,
            gamma_m = gamma_m,
            gamma_s = gamma_s,
            gamma_tilde_m = gmm_m,
            gamma_tilde_s = gmm_s,
            RI = ari_value,
            PosteriorModeClusters = post_mode_clusters,
            clusters_VI = length(cluster_ids)
          ))
        }
      }
    }
  }
}


# Save final results

file_name <- paste0("summary_simulation_results1_p",p,"_K",K,"_sigma",sigma_values[1],"_",sigma_values[2],"_",sigma_values[3],"totiter",totiter+burnin,"_n",sum(nj),".xlsx")


write_xlsx(results_summary, file_name)


cat("Simulation completed successfully!\n")


