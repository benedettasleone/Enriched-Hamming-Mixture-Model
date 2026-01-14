#include "labels.hpp"
#include <Rcpp.h>

void labels::compact_allocations(const int empty_index, bool inner, const int outer_index) {
    if (inner) {
        //find outer cluster index for the empty inner cluster
        if (outer_index == -1) {
        Rcpp::stop("outer_index must be provided when compacting inner allocations");
        }
        int outer_idx = outer_index;
        int max_inner = h_m[outer_idx] - 1;

        // riallocati ultimo inner cluster nell'empty cluster
        for(size_t i = 0; i < inner_allocations.size(); i++) 
            if (outer_allocations[i] == outer_idx && inner_allocations[i] == max_inner) 
                inner_allocations[i] = empty_index;

        // aggiorna h_m
        h_m[outer_idx]--;

        // aggiorno S
        S[outer_idx]--;

        // Sposto sigma + elimino ultimo
        inner_sigma[outer_idx][empty_index] = std::move(inner_sigma[outer_idx][max_inner]);
        inner_sigma[outer_idx].erase(max_inner);

        if (inner_weights.find(outer_idx) != inner_weights.end()) {
            if (inner_weights[outer_idx].find(max_inner) != inner_weights[outer_idx].end()) {
                inner_weights[outer_idx][empty_index] = inner_weights[outer_idx][max_inner];
                inner_weights[outer_idx].erase(max_inner);
            }
        }

    } else {
        // Trova il massimo outer
        int max_outer = -1;
        for (size_t i = 0; i < outer_allocations.size(); i++) {
            if (outer_allocations[i] > max_outer) {
                max_outer = outer_allocations[i];
            }
        }
        
        if (max_outer == -1 || max_outer == empty_index) {
            K--;
            M--;
            return;
        }

        // Rialloca
        for(size_t i = 0; i < outer_allocations.size(); i++) {
            if (outer_allocations[i] == max_outer) {
                outer_allocations[i] = empty_index;
            }
        }

        K--;
        M--;

        // Sposta outer_center
        outer_center[empty_index] = std::move(outer_center[max_outer]);
        outer_center.erase(max_outer);
        

        if (inner_sigma.find(max_outer) != inner_sigma.end()) {
            inner_sigma[empty_index] = std::move(inner_sigma[max_outer]);
            inner_sigma.erase(max_outer);  
        }
        
        // inner_weights
        if (inner_weights.find(max_outer) != inner_weights.end()) {
            inner_weights[empty_index] = std::move(inner_weights[max_outer]);
            inner_weights.erase(max_outer);  
        }
        
        // outer_weights
        if (outer_weights.find(max_outer) != outer_weights.end()) {
            outer_weights[empty_index] = outer_weights[max_outer];
            outer_weights.erase(max_outer);
        }
        
        if (max_outer < h_m.size() && empty_index < h_m.size()) {
            h_m[empty_index] = h_m[max_outer];
        }
        if (max_outer < S.size() && empty_index < S.size()) {
            S[empty_index] = S[max_outer];
        }
        
        // Rimuovi elementi oltre M
        h_m.resize(M);
        S.resize(M);
    }
}