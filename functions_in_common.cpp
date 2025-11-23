#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <random>
using namespace Rcpp;



double phi_dir(double u, double sig) {
  return std::pow(u + 1.0, -sig);
}


// Helper function: compute log factorial
double log_factorial(int n) {
  if (n <= 1) return 0.0;
  double result = 0.0;
  for (int i = 2; i <= n; i++) {
    result += std::log(static_cast<double>(i));
  }
  return result;
}


// Fast log factorial lookup
double fast_log_factorial(int x, NumericVector precomputed_log_factorial) {
  if (x < precomputed_log_factorial.size()) {
    return precomputed_log_factorial[x];
  }

  return lgamma(x + 1);
}


// [[Rcpp::export]]
int sampleM(int k, double u, IntegerVector x, double mu, double sig, NumericVector precomputed_log_factorial) {
  

  RNGScope rngScope;
  
  int n = x.size();
  NumericVector prob(n);
  
  // Pre-calculate common terms
  double log_phi_mu = std::log(phi_dir(u, sig) * mu);
  
  // Calculate log probabilities 
  double max_prob = -INFINITY;
  for (int i = 0; i < n; i++) {
    int xi = x[i];
    prob[i] = std::log(xi + k) - fast_log_factorial(xi, precomputed_log_factorial) + xi * log_phi_mu;
    if (prob[i] > max_prob) max_prob = prob[i];
  }
  
  // Convert to probabilities and normalize in one pass
  double sum_prob = 0.0;
  for (int i = 0; i < n; i++) {
    prob[i] = std::exp(prob[i] - max_prob);
    sum_prob += prob[i];
  }
  
  // Normalize
  for (int i = 0; i < n; i++) {
    prob[i] /= sum_prob;
  }
  
  // Sample from x using calculated probabilities
  int sampled_idx = Rcpp::sample(n, 1, false, prob)[0] - 1;
  return k + x[sampled_idx];
}


// [[Rcpp::export]]
IntegerVector sampleS(IntegerVector h_m, NumericVector u_low, IntegerVector x, 
                      double nu, double sig, NumericVector precomputed_log_factorial) {
  
  RNGScope rngScope;
  
  int M = h_m.size();
  int n = x.size();
  IntegerVector result(M);
  
  // Pre-compute log factorials
  NumericVector log_factorial_x(n);
  for (int i = 0; i < n; i++) {
    log_factorial_x[i] = fast_log_factorial(x[i], precomputed_log_factorial);
  }
  
  for (int mm = 0; mm < M; mm++) {
    
    if (h_m[mm] > 0) {
      NumericVector prob(n);
      
      // Pre-calculate terms for this iteration
      double log_phi_nu = std::log(phi_dir(u_low[mm], sig) * nu);
      int h_mm = h_m[mm];
      
      // Calculate log probabilities 
      double max_prob = -INFINITY;
      for (int i = 0; i < n; i++) {
        int xi = x[i];
        prob[i] = std::log(xi + h_mm) - log_factorial_x[i] + xi * log_phi_nu;
        if (prob[i] > max_prob) max_prob = prob[i];
      }
      
      // Convert to probabilities and normalize
      double sum_prob = 0.0;
      for (int i = 0; i < n; i++) {
        prob[i] = std::exp(prob[i] - max_prob);
        sum_prob += prob[i];
      }
      
      for (int i = 0; i < n; i++) {
        prob[i] /= sum_prob;
      }
      
      // Sample and compute result
      int sampled_idx = Rcpp::sample(n, 1, false, prob)[0] - 1;
      result[mm] = h_mm + x[sampled_idx];
      
    } else {
      // Sample from Poisson and add 1
      result[mm] = Rcpp::rpois(1, nu)[0] + 1;
    }
  }
  
  return result;
}



// [[Rcpp::export]]
NumericVector sample_w_un(int M, IntegerVector m, double u, double sig) {
  
  RNGScope rngScope;
  
  // Get unique values from m
  std::set<int> m_set(m.begin(), m.end());
  int unique_count = m_set.size();
  int empty = M - unique_count;
  
  // Count occurrences of each m value
  std::map<int, int> m_counts;
  for (int i = 0; i < m.size(); i++) {
    m_counts[m[i]]++;
  }
  
  // Create n vector: counts of existing clusters + zeros for empty clusters
  std::vector<double> n_vec;
  
  // Add counts for existing clusters
  std::vector<int> unique_m(m_set.begin(), m_set.end());
  std::sort(unique_m.begin(), unique_m.end());  
  
  for (int cluster : unique_m) {
    n_vec.push_back(static_cast<double>(m_counts[cluster]));
  }
  
  // Add zeros for empty clusters
  for (int i = 0; i < empty; i++) {
    n_vec.push_back(0.0);
  }
  
  // Sample from gamma distribution
  int n_length = n_vec.size();
  NumericVector w_un(n_length);
  
  for (int i = 0; i < n_length; i++) {
    double shape = n_vec[i] + sig;
    double rate = 1.0 + u;
    w_un[i] = Rcpp::rgamma(1, shape, 1.0/rate)[0];  
  }
  
  return w_un;
}



// [[Rcpp::export]]
NumericVector sample_q_un(int M, IntegerVector S, IntegerVector m, 
                         IntegerVector s, NumericVector u_low, double sig) {
  
  RNGScope rngScope;
  
  int total_S = std::accumulate(S.begin(), S.end(), 0);
  NumericVector q_un(total_S);
  
  for (int mm = 0; mm < M; mm++) {  
    
    // Calculate indices a and b 
    int a, b;
    if (mm == 0) {
      a = 0;
      b = S[0];
    } else {
      a = 0;
      for (int i = 0; i < mm; i++) {
        a += S[i];
      }
      b = a + S[mm];
    }
    
    // Get subset of s where m == (mm+1) 
    std::vector<int> s_subset;
    for (int i = 0; i < m.size(); i++) {
      if (m[i] == (mm + 1)) {  
        s_subset.push_back(s[i]);
      }
    }
    
    if (s_subset.size() > 0) {
      // Call sample_w_un
      IntegerVector s_subset_rcpp(s_subset.size());
      for (size_t i = 0; i < s_subset.size(); i++) {
        s_subset_rcpp[i] = s_subset[i];
      }
      
      NumericVector temp_w = sample_w_un(S[mm], s_subset_rcpp, u_low[mm], sig);
      
      for (int i = 0; i < temp_w.size() && i < (b - a); i++) {
        q_un[a + i] = temp_w[i];
      }
      
    } else {
      // Sample from gamma distribution
      for (int i = 0; i < (b - a); i++) {
        q_un[a + i] = Rcpp::rgamma(1, sig, 1.0)[0];
      }
    }
  }
  
  return q_un;
}


// [[Rcpp::export]]
NumericVector sample_u_low(int M, 
                           const IntegerVector& m, 
                           const NumericVector& q_un,
                           const IntegerVector& S) {
  
  // Create frequency table for m
  std::map<int, int> freq_table;
  for (int i = 0; i < m.size(); i++) {
    freq_table[m[i]]++;
  }
  
  // Get unique values and their counts
  std::vector<int> un_m;
  std::vector<int> n_m;
  for (const auto& pair : freq_table) {
    un_m.push_back(pair.first);
    n_m.push_back(pair.second);
  }
  
  int k = un_m.size();
  NumericVector sumq(k);  
  
  // Calculate sumq for each unique value
  for (int i = 0; i < k; i++) {
    int mm = un_m[i];
    int a, b;
    
    if (mm == 1) {
      a = 0;
      b = S[0]; 
    } else {
      
      a = 0;
      for (int j = 0; j < (mm - 1); j++) {
        a += S[j];
      }
      b = a + S[mm - 1]; 
    }
    
    
    sumq[i] = 0.0;
    for (int j = a; j < b; j++) {
      sumq[i] += q_un[j];
    }
  }
  
  
  NumericVector u_low(k);
  for (int i = 0; i < k; i++) {
    
    u_low[i] = R::rgamma(n_m[i], 1.0 / sumq[i]);
  }
  
  return u_low;
}




// [[Rcpp::export]]
List rearrange_ms(IntegerVector m, IntegerVector s, IntegerVector S, int M) {
  
  int n = m.size();
  IntegerVector m_temp = clone(m);
  IntegerVector S_temp = clone(S);
  IntegerVector s_temp = clone(s);
  
  // Get unique values from m
  std::set<int> m_set(m.begin(), m.end());
  std::vector<int> m_unique(m_set.begin(), m_set.end());
  int k = m_unique.size();
  
  IntegerVector h_m(k, 0);  // number of inner clusters per outer cluster
  
  // Count occurrences of each m value
  std::map<int, int> m_counts;
  for (int i = 0; i < n; i++) {
    m_counts[m[i]]++;
  }
  
  // Create old_m: sort by count (descending), then by value (ascending) for ties
  std::vector<int> old_m = m_unique;
  std::sort(old_m.begin(), old_m.end(), [&](int a, int b) {
    if (m_counts[a] != m_counts[b]) {
      return m_counts[a] > m_counts[b];  
    }
    return a < b;  
  });
  
  int count_m = 1;  
  int count_s = 1;  
  std::vector<int> stom;
  
  // Process each unique m value
  for (int x : old_m) {
    // Update m_temp
    for (int i = 0; i < n; i++) {
      if (m[i] == x) {
        m_temp[i] = count_m;
      }
    }
    
    S_temp[count_m - 1] = S[x - 1];  
    
    // Get unique s values for current m
    std::set<int> s_set_for_m;
    for (int i = 0; i < n; i++) {
      if (m[i] == x) {
        s_set_for_m.insert(s[i]);
      }
    }
    
    std::vector<int> s_unique(s_set_for_m.begin(), s_set_for_m.end());
    h_m[count_m - 1] = s_unique.size();
    
    // Count occurrences of each s value within current m
    std::map<int, int> s_counts;
    for (int i = 0; i < n; i++) {
      if (m[i] == x) {
        s_counts[s[i]]++;
      }
    }
    
    // Create old_s: sort by count (descending), then by value (ascending) for ties
    std::vector<int> old_s = s_unique;
    std::sort(old_s.begin(), old_s.end(), [&](int a, int b) {
      if (s_counts[a] != s_counts[b]) {
        return s_counts[a] > s_counts[b];  
      }
      return a < b;  
    });
    
    // Process each unique s value for current m
    for (int y : old_s) {
      // Update s_temp
      for (int i = 0; i < n; i++) {
        if (s[i] == y) {
          s_temp[i] = count_s;
        }
      }
      
      stom.push_back(count_m);
      count_s++;
    }
    
    // Handle case where S_temp[count_m] > h_m[count_m]
    if (S_temp[count_m - 1] > h_m[count_m - 1]) {
      int diff = S_temp[count_m - 1] - h_m[count_m - 1];
      for (int y = 0; y < diff; y++) {
        stom.push_back(count_m);
        count_s++;
      }
    }
    
    count_m++;
  }
  
  // Handle remaining m values not in m_unique (from 1 to M)
  std::set<int> processed_m(m_unique.begin(), m_unique.end());
  for (int x = 1; x <= M; x++) {
    if (processed_m.find(x) == processed_m.end()) {
      S_temp[count_m - 1] = S[x - 1];  
      
      for (int y = 0; y < S_temp[count_m - 1]; y++) {
        stom.push_back(count_m);
        count_s++;
      }
      
      count_m++;
    }
  }
  
  
  IntegerVector stom_rcpp(stom.size());
  for (size_t i = 0; i < stom.size(); i++) {
    stom_rcpp[i] = stom[i];
  }
  
  return List::create(
    Named("m") = m_temp,
    Named("s") = s_temp,
    Named("stom") = stom_rcpp,
    Named("S") = S_temp,
    Named("k") = k,
    Named("h_m") = h_m
  );
}

// [[Rcpp::export]]
List rearrange_ms2(IntegerVector m, IntegerVector s, IntegerVector S, int M, 
                   NumericMatrix center, IntegerVector attrsize, NumericMatrix zeta, double tao) {
  
  int n = m.size();
  int d = center.nrow();
  int total_S = std::accumulate(S.begin(), S.end(), 0);
  
  IntegerVector m_temp = clone(m);
  IntegerVector S_temp = clone(S);
  IntegerVector s_temp = clone(s);
  
  NumericMatrix center_new(d, total_S);
  std::fill(center_new.begin(), center_new.end(), 1.0);

  std::vector<bool> filled(total_S, false);
  
  std::set<int> m_set(m.begin(), m.end());
  std::vector<int> m_unique(m_set.begin(), m_set.end());
  int k = m_unique.size();
  
  IntegerVector h_m(k, 0);
  
  std::map<int, int> m_counts;
  for (int i = 0; i < n; i++) m_counts[m[i]]++;
  
  std::vector<int> old_m = m_unique;
  std::sort(old_m.begin(), old_m.end(), [&](int a, int b) {
    if (m_counts[a] != m_counts[b]) return m_counts[a] > m_counts[b];
    return a < b;
  });
  
  int count_m = 1;
  int count_s = 1;
  std::vector<int> stom;
  stom.reserve(total_S);
  
  for (int old_outer : old_m) {
    
    for (int i = 0; i < n; i++) if (m[i] == old_outer) m_temp[i] = count_m;
    
    if (count_m - 1 < S_temp.size()) {
      S_temp[count_m - 1] = S[old_outer - 1];
    } else {
      stop("rearrange_ms2: count_m-1 >= S_temp.size()");
    }
    
    std::set<int> s_set_for_m;
    for (int i = 0; i < n; i++) if (m[i] == old_outer) s_set_for_m.insert(s[i]);
    std::vector<int> s_unique(s_set_for_m.begin(), s_set_for_m.end());
    h_m[count_m - 1] = s_unique.size();
    
    std::map<int, int> s_counts;
    for (int i = 0; i < n; i++) if (m[i] == old_outer) s_counts[s[i]]++;
    
    std::vector<int> old_s = s_unique;
    std::sort(old_s.begin(), old_s.end(), [&](int a, int b) {
      if (s_counts[a] != s_counts[b]) return s_counts[a] > s_counts[b];
      return a < b;
    });
    
    for (int old_inner : old_s) {
      for (int i = 0; i < n; i++) if (s[i] == old_inner) s_temp[i] = count_s;
      
      stom.push_back(count_m);
      
      int old_col = old_inner - 1;
      int new_col = count_s - 1;
      if (new_col < 0 || new_col >= total_S) stop("rearrange_ms2: new_col fuori range (copy).");
      if (old_col >= 0 && old_col < center.ncol()) {
        for (int j = 0; j < d; j++) center_new(j, new_col) = center(j, old_col);
      }
      filled[new_col] = true;
      count_s++;
    }
    
    if (S_temp[count_m - 1] > h_m[count_m - 1]) {
      int diff = S_temp[count_m - 1] - h_m[count_m - 1];
      for (int y = 0; y < diff; y++) {
        stom.push_back(count_m);
        int new_col = count_s - 1;
        if (new_col < 0 || new_col >= total_S) stop("rearrange_ms2: new_col fuori range (empty cols).");
        if (!filled[new_col]) {
          for (int j = 0; j < d; j++) {
            int m_j = attrsize[j];
            if (m_j <= 0) stop("rearrange_ms2: attrsize must be > 0.");
            
            NumericVector prob(m_j);
            bool has_zeta = (count_m - 1) >= 0 && (count_m - 1) < zeta.ncol();
            
            if (has_zeta) {
              int zeta_val = (int) zeta(j, count_m - 1);
              for (int k2 = 0; k2 < m_j; k2++) {
                if ((k2 + 1) == zeta_val) prob[k2] = 1.0;
                else prob[k2] = std::exp(-1.0 / tao);
              }
            } else {
              for (int k2 = 0; k2 < m_j; k2++) {
                prob[k2] = 1.0 / m_j;
              }
            }
            
            double sum_prob = std::accumulate(prob.begin(), prob.end(), 0.0);
            if (sum_prob <= 0.0) {
              for (int k2 = 0; k2 < m_j; k2++) prob[k2] = 1.0 / m_j;
            } else {
              for (int k2 = 0; k2 < m_j; k2++) prob[k2] /= sum_prob;
            }
            
            IntegerVector levels = seq_len(m_j);
            center_new(j, new_col) = Rcpp::sample(levels, 1, false, prob)[0];
          }
          filled[new_col] = true;
        }
        count_s++;
      }
    }
    count_m++;
  }
  
  std::set<int> processed_m(m_unique.begin(), m_unique.end());
  for (int x = 1; x <= M; x++) {
    if (processed_m.find(x) == processed_m.end()) {
      if (count_m - 1 < 0 || count_m - 1 >= S_temp.size()) stop("rearrange_ms2: count_m-1 fuori range in unobserved block.");
      S_temp[count_m - 1] = S[x - 1];
      
      for (int y = 0; y < S_temp[count_m - 1]; y++) {
        stom.push_back(count_m);
        int new_col = count_s - 1;
        if (new_col < 0 || new_col >= total_S) stop("rearrange_ms2: new_col fuori range (unobserved outer).");
        if (!filled[new_col]) {
          for (int j = 0; j < d; j++) {
            int m_j = attrsize[j];
            if (m_j <= 0) stop("rearrange_ms2: attrsize must be > 0.");
            
            NumericVector prob(m_j);
            bool has_zeta = (count_m - 1) >= 0 && (count_m - 1) < zeta.ncol();
            
            if (has_zeta) {
              int zeta_val = (int) zeta(j, count_m - 1);
              for (int k2 = 0; k2 < m_j; k2++) {
                if ((k2 + 1) == zeta_val) prob[k2] = 1.0;
                else prob[k2] = std::exp(-1.0 / tao);
              }
            } else {
              for (int k2 = 0; k2 < m_j; k2++) {
                prob[k2] = 1.0 / m_j;
              }
            }
            
            double sum_prob = std::accumulate(prob.begin(), prob.end(), 0.0);
            if (sum_prob <= 0.0) {
              for (int k2 = 0; k2 < m_j; k2++) prob[k2] = 1.0 / m_j;
            } else {
              for (int k2 = 0; k2 < m_j; k2++) prob[k2] /= sum_prob;
            }
            
            IntegerVector levels = seq_len(m_j);
            center_new(j, new_col) = Rcpp::sample(levels, 1, false, prob)[0];
          }
          filled[new_col] = true;
        }
        count_s++;
      }
      count_m++;
    }
  }
  
  if ((int)stom.size() > total_S) {
    Rcpp::Rcout << "Warning: stom.size() (" << stom.size() << ") > total_S (" << total_S << "). Some columns ignored." << std::endl;
  }
  
  IntegerVector stom_rcpp(stom.size());
  for (size_t i = 0; i < stom.size(); i++) stom_rcpp[i] = stom[i];
  
  return List::create(
    Named("m") = m_temp,
    Named("s") = s_temp,
    Named("stom") = stom_rcpp,
    Named("S") = S_temp,
    Named("k") = k,
    Named("h_m") = h_m,
    Named("center") = center_new
  );
}

