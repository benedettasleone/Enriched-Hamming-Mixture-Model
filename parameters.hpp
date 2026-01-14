#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

struct parameters {
    NumericVector v;
    NumericVector w;
    double sig_m;
    double sig_s;
    double mu;
    double nu;
    int d; // number of variables
    IntegerVector attrsize; // sizes of categorical variables
    NumericMatrix data; // data matrix

    std::vector<double> precomputed_log_factorial; // Precomputed log factorials for efficiency
    
    int n_probs_poisson = 1000; // Predefined size for Poisson probabilities

    // Default constructor
    parameters() : sig_m(0.0), sig_s(0.0), mu(0.0), nu(0.0), d(0), n_probs_poisson(1000) {}
    
    // Parameterized constructor
    parameters(NumericVector v_, NumericVector w_, double sig_m_, double sig_s_, 
               double mu_, double nu_, int d_, IntegerVector attrsize_, NumericMatrix data_)
        : v(v_), w(w_), sig_m(sig_m_), sig_s(sig_s_), mu(mu_), nu(nu_), 
          d(d_), attrsize(attrsize_), data(data_), n_probs_poisson(1000) {
        
        // Precompute log factorials for efficiency
        precomputed_log_factorial.resize(n_probs_poisson);
        precomputed_log_factorial[0] = 0.0;
        for (int i = 1; i < n_probs_poisson; i++) {
            precomputed_log_factorial[i] = precomputed_log_factorial[i-1] + std::log(i);
        }
    }

    
};


#endif
