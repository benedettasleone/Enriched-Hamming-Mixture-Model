#pragma once

#include <Rcpp.h>
#include "parameters.hpp"
#include "labels.hpp"
using namespace Rcpp;

class Sampler {
    private:
        const parameters& params;
        labels& lbls;
        double u_up; // variabile ausiliaria per il campionamento di M
        NumericVector u_low; // variabile ausiliaria per il campionamento di S

        RNGScope rngScope; // Initialize RNG once for the Sampler instance

        void sample_inner_allocations();
        void sample_M();
        void sample_S();
        void sample_centers();
        void sample_sigmas();
        void sample_inner_weights();
        void sample_outer_weights();
        void sample_u_up();
        void sample_u_low();

        // Helper functions
        double phi_dir(double sig);
        double fast_log_factorial(int x, const std::vector<double>& precomputed_log_factorial);
        std::unordered_map<int, double> sample_weights(int ncls, const std::vector<int>& alloc, double u, double sig);
        NumericMatrix kernel_hamming_rcpp_matrix(NumericMatrix y, std::unordered_map<int, std::vector<double>> center, 
            std::unordered_map<int, std::unordered_map<int, std::vector<double>>> sigma, IntegerVector attrsize, int d, int total_S, 
            std::map<int, std::pair<int, int>> inner_allocation_mapping);

        
    public: 
        Sampler(const parameters& parameters, labels& labels)
            : params(parameters), lbls(labels) {}

        void sample();


};