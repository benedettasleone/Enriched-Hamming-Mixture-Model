#include <Rcpp.h>
#include "sampler.hpp"
#include "hyperg.hpp"
#include "parameters.hpp"
#include "labels.hpp"
#include <chrono>
#include <unordered_map>
#include <vector>
using namespace Rcpp;


// [[Rcpp::export]]
List run_mcmc_cpp(NumericMatrix Y, IntegerVector m_init, IntegerVector s_init,
                  NumericVector w_un_init, NumericVector q_un_init,
                  NumericMatrix center_init, NumericMatrix sigma_init,
                  int M_init, IntegerVector S_init, IntegerVector attrsize, NumericVector v_init,
                  NumericVector w_init, int burnin, int totiter, double mu, double nu, 
                  double sig_m, double sig_s, bool verbose = false) {
    
    // Initialize RNG
    RNGScope rngScope;
    
    // Input validation
    int n = Y.nrow();
    int p = Y.ncol();
    
    if (m_init.size() != n || s_init.size() != n) {
        stop("Dimension mismatch: m_init and s_init must have length n");
    }
    if (center_init.nrow() != p) {
        stop("center_init must have p rows");
    }
    if (sigma_init.nrow() != p) {
        stop("sigma_init must have p rows");
    }
    
    if (verbose) {
        Rcout << "Starting MCMC with n=" << n << ", p=" << p << std::endl;
        Rcout << "burnin=" << burnin << ", totiter=" << totiter << std::endl;
    }
    
    // Initialize variables
    std::vector<int> m(m_init.begin(), m_init.end());
    std::vector<int> s(s_init.begin(), s_init.end());
    // from w_un_init to un unordered_map
    std::unordered_map<int, double> w_un;
    for(int i = 0; i < w_un_init.size(); i++) {
        w_un[i] = w_un_init[i];
    }
    // from q_un_init to un unordered_map of unordered_map
    std::unordered_map<int, std::unordered_map<int, double>> q_un;
    int current_inner_idx = 0;
    for(int mm = 0; mm < M_init; mm++) {
        int num_inner_clusters = S_init[mm];
        for(int ss = 0; ss < num_inner_clusters; ss++) {
            q_un[mm][ss] = q_un_init[current_inner_idx];
            current_inner_idx++;
        }
    }
    //from sigma_init to unordered_map, not all outer have same number of inner clusters, so we need to build the map accordingly
    std::unordered_map<int, std::unordered_map<int, std::vector<double>>> inner_sigma;
    current_inner_idx = 0;
    for(int mm = 0; mm < M_init; mm++) {
        int num_inner_clusters = S_init[mm];
        for(int ss = 0; ss < num_inner_clusters; ss++) {
            std::vector<double> sigma_vec;
            for(int j = 0; j < p; j++) {
                sigma_vec.push_back(sigma_init(j, current_inner_idx));
            }
            inner_sigma[mm][ss] = sigma_vec;
            current_inner_idx++;
        }
    }
    //from center_init to unordered_map
    std::unordered_map<int, std::vector<double>> outer_center;
    for(int mm = 0; mm < M_init; mm++) {
        std::vector<double> center_vec;
        for(int j = 0; j < p; j++) {
            center_vec.push_back(center_init(j, mm));
        }
        outer_center[mm] = center_vec;
    }
    int M = M_init;
    std::vector<int> S(S_init.begin(), S_init.end());

    
    // Create parameters and labels objects
    parameters params(v_init, w_init, sig_m, sig_s, mu, nu, p, attrsize, Y);
    labels lbls(s, m, inner_sigma, outer_center, w_un, q_un, M, S);

    // Pre-allocate storage for post-burnin samples 
    NumericMatrix m_results(totiter, n);
    NumericMatrix s_results(totiter, n);
    List center_results(totiter);
    List sigma_results(totiter);
    IntegerVector M_results(totiter);
    IntegerVector k_results(totiter);
    List S_results(totiter);
    List h_m_results(totiter);
    List w_un_results(totiter);
    List q_un_results(totiter);
    NumericVector u_up_results(totiter);
    List u_low_results(totiter);
   
    
    // Main MCMC loop
    int post_burnin_idx = 0;
    
    for (int iter = 0; iter < (burnin + totiter); iter++) {
        
       
        if (verbose && iter == 0) {
            Rcout << "MCMC has started" << std::endl;
        }
        if (verbose && (iter + 1) % 500 == 0) {  
            Rcout << "Iteration: " << (iter + 1) << std::endl;
        }

        auto start = std::chrono::high_resolution_clock::now();
        if(iter == 0 or iter == (burnin + totiter -1) or iter == burnin){
          start = std::chrono::high_resolution_clock::now();
        }

            Sampler sampler(params, lbls);
            sampler.sample();
            
            
            // Store post-burnin samples 
            if (iter >= burnin) {
            
            for (int i = 0; i < n; i++) {
                m_results(post_burnin_idx, i) = lbls.get_outer_allocations()[i];
                s_results(post_burnin_idx, i) = lbls.get_inner_allocations()[i];
            }
            
            // Convert center back to matrix
            NumericMatrix center_mat(p, lbls.get_M());
            for(int mm = 0; mm < lbls.get_M(); mm++) {
                for(int j = 0; j < p; j++) {
                    center_mat(j, mm) = lbls.get_outer_center().at(mm)[j];
                }
            }
            center_results[post_burnin_idx] = center_mat;
            
            // Convert sigma back to matrix
            int total_S = std::accumulate(lbls.get_S().begin(), lbls.get_S().end(), 0);
            NumericMatrix sigma_mat(p, total_S);
            int global_ss = 0;
            for(int mm = 0; mm < lbls.get_M(); mm++) {
                for(int ss = 0; ss < lbls.get_S()[mm]; ss++) {
                    for(int j = 0; j < p; j++) {
                        sigma_mat(j, global_ss) = lbls.get_inner_sigma().at(mm).at(ss)[j];
                    }
                    global_ss++;
                }
            }
            sigma_results[post_burnin_idx] = sigma_mat;
            
            M_results[post_burnin_idx] = lbls.get_M();
            k_results[post_burnin_idx] = lbls.get_K();
            S_results[post_burnin_idx] = wrap(lbls.get_S());
            h_m_results[post_burnin_idx] = wrap(lbls.get_h_m());
            w_un_results[post_burnin_idx] = wrap(lbls.get_outer_weights());
            q_un_results[post_burnin_idx] = wrap(lbls.get_inner_weights());
            
            post_burnin_idx++;
        }
        
        if(iter == 0 or iter == (burnin + totiter -1) or iter == burnin){
          auto end = std::chrono::high_resolution_clock::now(); 
          std::chrono::duration<double> elapsed = end - start;
          Rcout << "Time for iteration " << (iter + 1) << ": " << elapsed.count() << " seconds" << std::endl;
        }
            
        
        
    }
    
    if (verbose) {
        Rcout << "MCMC completed successfully!" << std::endl;
    }
    
    // Build final result list once at the end
    List post_burnin_chain = List::create(
        Named("m") = m_results,
        Named("s") = s_results,
        Named("center") = center_results,
        Named("sigma") = sigma_results,
        Named("M") = M_results, 
        Named("k") = k_results,
        Named("S") = S_results,
        Named("S_m") = h_m_results,
        Named("w_un") = w_un_results,
        Named("q_un") = q_un_results,
        Named("u_up") = u_up_results,
        Named("u_low") = u_low_results
    );
    
    return post_burnin_chain;
}





