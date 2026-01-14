#include "sampler.hpp"
#include "hyperg.hpp"


NumericMatrix Sampler::kernel_hamming_rcpp_matrix(NumericMatrix y, std::unordered_map<int, std::vector<double>> center, std::unordered_map<int, std::unordered_map<int, std::vector<double>>> sigma, IntegerVector attrsize, int d, int total_S, std::map<int, std::pair<int, int>> inner_allocation_mapping) {
  int n = y.nrow();
  NumericMatrix out(n, total_S);

  for (int i = 0; i < n; i++){
    for(int ss = 0; ss < total_S; ss++){
      std::pair<int, int> outer_inner_indices = inner_allocation_mapping[ss];
      int mm = outer_inner_indices.first; 
      int inner_cluster_local_idx = outer_inner_indices.second;

      double val = 0.0;
      for(int j = 0; j < d; j++) {
        int diff = (y(i, j) == center[mm][j]) ? 0 : 1;
        double numerator = - diff / sigma[mm][inner_cluster_local_idx][j];
        double exp_term = std::exp(1.0 / sigma[mm][inner_cluster_local_idx][j]);
        double attr_ratio = (attrsize[j] - 1.0) / exp_term;
        double denominator = std::log(1.0 + attr_ratio);
        val += numerator - denominator;
      }
      out(i, ss) = val;
    }
  }
  return out;
}

void Sampler::sample_inner_allocations() {

  std::vector<int> S = lbls.get_S();
  int total_S = std::accumulate(S.begin(), S.end(), 0);
  
  //now inner allocation count starts from 0 for each outer cluster, so we need to map it back. I do it here so it's done only once. inner allocation_mapping should tell me for each global inner cluster index which outer cluster and inner cluster is, cause I need for each inner cluster its outer cluster and it's index within that outer cluster
  std::map<int, std::pair<int, int>> inner_allocation_mapping; // key: global inner cluster index, value: (outer cluster index, inner cluster index within outer cluster)
  int inner_cluster_idx = 0;
  for (int mm = 0; mm < lbls.get_M(); mm++) {
    for (int ss = 0; ss < S[mm]; ss++) {
      inner_allocation_mapping[inner_cluster_idx] = std::make_pair(mm, ss);
      inner_cluster_idx++;
    }
  }

  // Compute kernel matrix: n x total_S
  NumericMatrix kernel_mat = kernel_hamming_rcpp_matrix(params.data, lbls.get_outer_center(), lbls.get_inner_sigma(), params.attrsize, params.d, total_S, inner_allocation_mapping);

  // Compute log weights vector for each inner cluster
  NumericVector log_wq(total_S);
  int start_idx = 0; 
  
  for (int mm = 0; mm < lbls.get_M(); mm++) {
    double log_w = std::log(lbls.get_outer_weights().at(mm));
    
    for (int ss = 0; ss < S[mm]; ss++) {
      int global_idx = start_idx + ss;
      log_wq[global_idx] = log_w + std::log(lbls.get_inner_weights().at(mm).at(ss));
    }
    start_idx += S[mm];
  }

  for (int i = 0; i < params.data.nrow(); i++) {
    // Vectorized log probabilities for all inner clusters for observation i
    NumericVector log_prob(total_S);
    for (int j = 0; j < total_S; j++) {
      log_prob[j] = kernel_mat(i, j) + log_wq[j];
    }
    
    // Stabilize log probabilities to avoid numerical issues
    NumericVector exp_log_prob(total_S);
    for (int j = 0; j < total_S; j++) {
      exp_log_prob[j] = std::exp(log_prob[j]);
    }
    
    double sum_exp = std::accumulate(exp_log_prob.begin(), exp_log_prob.end(), 0.0);
    
    if (sum_exp == 0.0) {
    double min_log_prob = *std::min_element(log_prob.begin(), log_prob.end());
    for (int j = 0; j < total_S; j++) {
      log_prob[j] -= min_log_prob;
    }
    } else if (std::isinf(sum_exp)) {
      double max_log_prob = *std::max_element(log_prob.begin(), log_prob.end());
      for (int j = 0; j < total_S; j++) {
        log_prob[j] -= max_log_prob;
      }
    }

    // Convert to probabilities
    NumericVector final_probs(total_S);
    for (int j = 0; j < total_S; j++) {
      final_probs[j] = std::exp(log_prob[j]);
    }
    
    double prob_sum = std::accumulate(final_probs.begin(), final_probs.end(), 0.0);
    for (int j = 0; j < total_S; j++) {
      final_probs[j] /= prob_sum;
    }

    // Sample s[i] from inner clusters weighted by final_probs
    IntegerVector indices = seq_len(total_S);
    int sampled_idx = Rcpp::sample(indices, 1, false, final_probs)[0] - 1;

    //now inner allocation count starts from 0 for each outer cluster, so we need to map it back
    std::pair<int, int> outer_inner_indices = inner_allocation_mapping[sampled_idx];
    int outer_cluster_idx = outer_inner_indices.first;
    int inner_cluster_idx = outer_inner_indices.second;

    lbls.set_outer_allocation_at(i, outer_cluster_idx);
    lbls.set_inner_allocation_at(i, inner_cluster_idx);
  }
  
}


double Sampler::phi_dir(double sig) {
  return std::pow(this->u_up + 1.0, -sig);
}


double Sampler::fast_log_factorial(int x, const std::vector<double>& precomputed_log_factorial) {
  if (x < precomputed_log_factorial.size()) {
    return precomputed_log_factorial[x];
  }

  return lgamma(x + 1);
}


void Sampler::sample_M() {
  
  int n_prb = params.n_probs_poisson;
  NumericVector prob(n_prb);

  int old_K = lbls.get_K();
  
  // Pre-calculate common terms
  double log_phi_mu = std::log(phi_dir(params.sig_m) * params.mu);
  
  // Calculate log probabilities with optimizations
  double max_prob = -INFINITY;
  for (int i = 0; i < n_prb; i++) {
    int xi = i;
    prob[i] = std::log(xi + old_K) - fast_log_factorial(xi, params.precomputed_log_factorial) + xi * log_phi_mu;
    if (prob[i] > max_prob) max_prob = prob[i];
  }
  
  // Convert to probabilities and normalize in one pass
  double sum_prob = 0.0;
  for (int i = 0; i < n_prb; i++) {
    prob[i] = std::exp(prob[i] - max_prob);
    sum_prob += prob[i];
  }
  
  // Normalize
  for (int i = 0; i < n_prb; i++) {
    prob[i] /= sum_prob;
  }
  
  // Sample from x using calculated probabilities
  int sampled_idx = Rcpp::sample(n_prb, 1, false, prob)[0] - 1;
  lbls.set_M( old_K + sampled_idx );

}


void Sampler::sample_S() {
    
  int n_prb = params.n_probs_poisson;
  std::vector<int> result(lbls.get_M());

  // Pre-calculate common terms
  NumericVector log_factorial_x(n_prb);
  for (int i = 0; i < n_prb; i++) {
    log_factorial_x[i] = fast_log_factorial(i, params.precomputed_log_factorial);
  }

  //h_m extended
  IntegerVector h_m_extended(lbls.get_M());
  for (int i = 0; i < lbls.get_M(); i++) {
    h_m_extended[i] = lbls.get_h_m()[i];
  }
  for (int i = lbls.get_h_m().size(); i < lbls.get_M(); i++) {
    h_m_extended[i] = 0;
  }
  
  for (int mm = 0; mm < lbls.get_M(); mm++) {

    if (h_m_extended[mm] > 0) {
      NumericVector prob(n_prb);

      // Pre-calculate terms for this iteration
      double log_phi_nu = std::log(phi_dir(params.sig_s) * params.nu);
      int h_mm = h_m_extended[mm];
      
      // Calculate log probabilities 
      double max_prob = -INFINITY;
      for (int i = 0; i < n_prb; i++) {
        int xi = i;
        prob[i] = std::log(xi + h_mm) - log_factorial_x[i] + xi * log_phi_nu;
        if (prob[i] > max_prob) max_prob = prob[i];
      }
      
      // Convert to probabilities and normalize
      double sum_prob = 0.0;
      for (int i = 0; i < n_prb; i++) {
        prob[i] = std::exp(prob[i] - max_prob);
        sum_prob += prob[i];
      }
      
      for (int i = 0; i < n_prb; i++) {
        prob[i] /= sum_prob;
      }
      
      // Sample and compute result
      int sampled_idx = Rcpp::sample(n_prb, 1, false, prob)[0] - 1;
      result[mm] = h_mm + sampled_idx;
    }
    else {
      // Sample from Poisson and add 1
      result[mm] = Rcpp::rpois(1, params.nu)[0] + 1;
    }
  }

  lbls.set_S(result);  
}


std::unordered_map<int, double> Sampler::sample_weights(int ncls, const std::vector<int>& alloc, double u, double sig) {
  
  //get unique clusters
  std::set<int> alloc_set(alloc.begin(), alloc.end());
  int unique_count = alloc_set.size();
  int empty = ncls - unique_count;
  
  // Count occurrences of each allocation value 
  std::map<int, int> alloc_counts;
  for (size_t i = 0; i < alloc.size(); i++) {
    alloc_counts[alloc[i]]++;
  }

  // Create n vector: counts of existing clusters + zeros for empty clusters
  std::vector<double> n_vec;
  
  // Add counts for existing clusters 
  std::vector<int> unique_alloc(alloc_set.begin(), alloc_set.end());
  std::sort(unique_alloc.begin(), unique_alloc.end());  

  for (int cluster : unique_alloc) {
    n_vec.push_back(static_cast<double>(alloc_counts[cluster]));
  }
  
  // Add zeros for empty clusters
  for (int i = 0; i < empty; i++) {
    n_vec.push_back(0.0);
  }
  
  // Sample from gamma distribution: rgamma(length(n), n + sig, 1 + u)
  int n_length = n_vec.size();
  std::unordered_map<int, double> weights;

  for (int i = 0; i < n_length; i++) {
    double shape = n_vec[i] + sig;
    double rate = 1.0 + u;
    weights[i] = Rcpp::rgamma(1, shape, 1.0/rate)[0];  
  }

  return weights;
}


void Sampler::sample_outer_weights(){
  lbls.set_outer_weights(sample_weights(lbls.get_M(), lbls.get_outer_allocations(), u_up, params.sig_m));
}


void Sampler::sample_inner_weights(){

  std::vector<int> S = lbls.get_S();

  std::unordered_map<int, std::unordered_map<int, double>> q_un;
  
  for (int mm = 0; mm < lbls.get_M(); mm++) { 
    
    // Get subset of s where m == (mm+1) 
    std::vector<int> s_subset;
    for (int i = 0; i < params.data.nrow(); i++) {
      if (lbls.get_outer_allocations()[i] == (mm)) {  
        s_subset.push_back(lbls.get_inner_allocations()[i]);
      }
    }

    if (s_subset.size() > 0) {
      // Do sample_inner_weights with the subset
      std::vector<int> s_subset_rcpp(s_subset.size());
      for (size_t i = 0; i < s_subset.size(); i++) {
        s_subset_rcpp[i] = s_subset[i];
      }

      std::unordered_map<int, double> temp_w = sample_weights(S[mm], s_subset_rcpp, u_low[mm], params.sig_s);

      for (const auto& [ss, weight] : temp_w) {
        q_un[mm][ss] = weight;
      }

      
    } else {
      // Sample from gamma distribution: rgamma(S[mm], sig, 1)
      
      for (int i = 0; i < S[mm]; i++) {
        q_un[mm][i] = Rcpp::rgamma(1, params.sig_s, 1.0)[0];
      }
    } 
  }

  lbls.set_inner_weights(q_un);

}


void Sampler::sample_u_up() {

  double w_un_sum = 0.0;
  for (const auto& [mm, weight] : lbls.get_outer_weights()) {
    w_un_sum += weight;
  }
  u_up = R::rgamma(params.data.nrow(), 1.0 / w_un_sum);

}


void Sampler::sample_u_low() {

  u_low = NumericVector(lbls.get_M());

  // Create frequency table for m
  std::map<int, int> freq_table;
  for (size_t i = 0; i < lbls.get_outer_allocations().size(); i++) {
    int mm = lbls.get_outer_allocations()[i];
    freq_table[mm]++;
  }
  
  std::vector<int> S = lbls.get_S();
  
  // Calculate sumq for each unique value
  for (const auto& [mm, n_mm] : freq_table) {
    double sum_inner_weights = 0.0;
    for (int local_ss = 0; local_ss < S[mm]; local_ss++) {
      sum_inner_weights += lbls.get_inner_weights().at(mm).at(local_ss);
    }
    
    // Sample u_low[mm] from Gamma(n_mm, 1/sum_inner_weights)
    u_low[mm] = R::rgamma(n_mm, 1.0 / sum_inner_weights);
  }
}


void Sampler::sample_centers() {

  int M = lbls.get_M();
  std::vector<int> h_m = lbls.get_h_m();

  std::unordered_map<int, std::vector<double>> c; // key: outer cluster mm, value: center vector for that cluster

  for (int mm = 0; mm < M; mm++) { 

    c[mm] = std::vector<double>();
    
    for (int j = 0; j < params.d; j++) {   
      
      int m_j = params.attrsize[j];
      std::vector<double> freq(m_j, 0.0);
      std::vector<double> prob(m_j, 0.0);
      
      for (int ms = 0; ms < h_m[mm]; ms++) {  
        
        // Count observations and frequencies for attribute j in cluster ms
        int n_ms = 0;
        std::fill(freq.begin(), freq.end(), 0.0);
        for (int i = 0; i < params.data.nrow(); i++) {
          if (lbls.get_outer_allocations()[i] == mm && lbls.get_inner_allocations()[i] == ms) {  
            n_ms++;
            int y_val = params.data(i, j) - 1;  
            if (y_val >= 0 && y_val < m_j) {
              freq[y_val] += 1.0;
            }
          }
        }

        
        // Update probabilities
        for (int k = 0; k < m_j; k++) {
          prob[k] -= (n_ms - freq[k]) / lbls.get_inner_sigma().at(mm).at(ms)[j];
        }
        
      }


      // Numerical stability: subtract max and exp
      double max_val = *std::max_element(prob.begin(), prob.end());
      for (int k = 0; k < m_j; k++) {
        prob[k] = std::exp(prob[k] - max_val);
      }

      // Normalize probabilities
      double sum_prob = std::accumulate(prob.begin(), prob.end(), 0.0);
      for (int k = 0; k < m_j; k++) {
        prob[k] /= sum_prob;
      }
      
      // Sample from 1:m_j with calculated probabilities
      NumericVector prob_rcpp(m_j);
      for (int k = 0; k < m_j; k++) {
        prob_rcpp[k] = prob[k];
      }

      IntegerVector choices = seq_len(m_j);
      int sampled_val = Rcpp::sample(choices, 1, false, prob_rcpp)[0];
      c[mm].push_back(sampled_val);

    }
  }

  lbls.set_outer_center(c);

}


void Sampler::sample_sigmas() {

  std::vector<int> S = lbls.get_S();
  int total_S = std::accumulate(S.begin(), S.end(), 0);
  
  std::unordered_map<int, std::unordered_map<int, std::vector<double>>> inner_sigma; // chiave esterna: outer cluster, chiave interna: inner cluster locale, valore: vettore di sigma per ogni dimensione
  NumericMatrix new_v(params.d, total_S);
  NumericMatrix new_w(params.d, total_S);

  // map each global inner cluster ss to its outer cluster mm and local inner cluster index
  std::map<int, std::pair<int, int>> inner_allocation_mapping;
  int inner_cluster_idx = 0;
  for (int mm = 0; mm < lbls.get_M(); mm++) {
    for (int ss = 0; ss < S[mm]; ss++) {
      inner_allocation_mapping[inner_cluster_idx] = std::make_pair(mm, ss);
      inner_cluster_idx++;
    }
  }

  for (int ss = 0; ss < total_S; ss++) {  

    std::vector<int> idx_ss;
    for (int i = 0; i < params.data.nrow(); i++) {
      if (lbls.get_inner_allocations()[i] == inner_allocation_mapping[ss].second && lbls.get_outer_allocations()[i] == inner_allocation_mapping[ss].first) {  
        idx_ss.push_back(i);
      }
    }
    
    int n = idx_ss.size();
    int mm = inner_allocation_mapping[ss].first;

    std::vector<int> match_counts(params.d, 0);

    // Count matches for each attribute
    for (int i : idx_ss) {
      for (int j = 0; j < params.d; j++) {
        if (params.data(i, j) == lbls.get_outer_center().at(mm)[j]) {
          match_counts[j]++;
        }
      }
    }

    for (int j = 0; j < params.d; j++) {
      int m_j = params.attrsize[j];
      int sumdelta = match_counts[j];

      new_w(j, ss) = params.w[j] + n - sumdelta;
      new_v(j, ss) = params.v[j] + sumdelta;
      
      // Sample from the full conditional
      double threshold = static_cast<double>(m_j - 1) / m_j;
      double qbeta_val = R::qbeta(0.1, new_w(j, ss) + 1, new_v(j, ss) - 1, 1, 0);
      
      double out;
      if (qbeta_val < threshold && threshold > 0.8) {
        // Use rejection sampling with beta distribution
        double x;
        do {
          x = R::rbeta(new_w(j, ss) + 1, new_v(j, ss) - 1);
        } while (x > threshold);
        
        out = x / ((m_j - 1) * (1 - x));
      } else {
        // Use bisection method
        double Omega = R::runif(0.0, 1.0);
        out = bisec_hyper2(new_w(j, ss), new_v(j, ss), m_j, Omega);
      }
      
      int local_ss = inner_allocation_mapping[ss].second;

      if (inner_sigma.find(mm) == inner_sigma.end()) {
        inner_sigma[mm] = std::unordered_map<int, std::vector<double>>();
      }
      if (inner_sigma[mm].find(local_ss) == inner_sigma[mm].end()) {
        inner_sigma[mm][local_ss] = std::vector<double>();
      }

      inner_sigma[mm][local_ss].push_back(-1.0 / std::log(out));

    }
  }

  lbls.set_inner_sigma(inner_sigma);

}


void Sampler::sample() {
  sample_inner_allocations();
  sample_u_up();
  sample_M();
  sample_outer_weights();
  sample_centers();
  sample_u_low();
  sample_S();
  sample_inner_weights();
  sample_sigmas();
}