#include "labels.hpp"
#include <Rcpp.h>

void labels::compact_allocations(const int empty_index, bool inner, const int outer_index) {
    
    if (inner) {
        // INNER clusters compaction
        
        if (outer_index == -1) {
            Rcpp::stop("outer_index must be provided when compacting inner allocations");
        }
        int outer_idx = outer_index;
        
        // Find the TRUE maximum inner index for this outer
        int max_inner = -1;
        for (size_t i = 0; i < inner_allocations.size(); i++) {
            if (outer_allocations[i] == outer_idx) {
                if (inner_allocations[i] > max_inner) {
                    max_inner = inner_allocations[i];
                }
            }
        }
        
        // Impossible case: no cluster allocated for this outer
        if (max_inner == -1) {
            h_m[outer_idx] = 0;
            return;
        }

        // If empty_index is already the last one, only remove the structures
        if (max_inner == empty_index) {
            h_m[outer_idx]--;
            inner_sigma[outer_idx].erase(empty_index);
            if (inner_weights.find(outer_idx) != inner_weights.end()) {
                inner_weights[outer_idx].erase(empty_index);
            }
            return;
        }

        // Normal case: move max_inner → empty_index
        
        // Reallocate observations
        for(size_t i = 0; i < inner_allocations.size(); i++) {
            if (outer_allocations[i] == outer_idx && inner_allocations[i] == max_inner) {
                inner_allocations[i] = empty_index;
            }
        }

        // Update h_m
        h_m[outer_idx]--;

        // Move sigma + delete last
        if (inner_sigma[outer_idx].find(max_inner) != inner_sigma[outer_idx].end()) {
            inner_sigma[outer_idx][empty_index] = std::move(inner_sigma[outer_idx][max_inner]);
            inner_sigma[outer_idx].erase(max_inner);
        }

        if (inner_weights.find(outer_idx) != inner_weights.end()) {
            if (inner_weights[outer_idx].find(max_inner) != inner_weights[outer_idx].end()) {
                inner_weights[outer_idx][empty_index] = inner_weights[outer_idx][max_inner];
                inner_weights[outer_idx].erase(max_inner);
            }
        }

    } else {
        // OUTER clusters compaction
        
        // Find the maximum outer
        int max_outer = -1;
        for (size_t i = 0; i < outer_allocations.size(); i++) {
            if (outer_allocations[i] > max_outer) {
                max_outer = outer_allocations[i];
            }
        }
        
        // Impossible case: no outer allocated
        if (max_outer == -1) {
            K = 0;
            return;
        }
        
        // If empty_index is already the last one, only remove the structures
        if (max_outer == empty_index) {
            K--;
            
            // Remove all structures for empty_index
            outer_center.erase(empty_index);
            outer_weights.erase(empty_index);
            inner_sigma.erase(empty_index);
            inner_weights.erase(empty_index);
            
            // Zero only h_m (S remains unchanged, it's a model parameter)
            if (empty_index < h_m.size()) {
                h_m[empty_index] = 0;
            }
            
            return;
        }

        // Normal case: move max_outer → empty_index
        
        // Reallocate observations
        for(size_t i = 0; i < outer_allocations.size(); i++) {
            if (outer_allocations[i] == max_outer) {
                outer_allocations[i] = empty_index;
            }
        }

        K--;

        // Move all structures
        if (outer_center.find(max_outer) != outer_center.end()) {
            outer_center[empty_index] = std::move(outer_center[max_outer]);
            outer_center.erase(max_outer);
        }

        if (inner_sigma.find(max_outer) != inner_sigma.end()) {
            inner_sigma[empty_index] = std::move(inner_sigma[max_outer]);
            inner_sigma.erase(max_outer);
        }
        
        if (inner_weights.find(max_outer) != inner_weights.end()) {
            inner_weights[empty_index] = std::move(inner_weights[max_outer]);
            inner_weights.erase(max_outer);
        }
        
        if (outer_weights.find(max_outer) != outer_weights.end()) {
            outer_weights[empty_index] = outer_weights[max_outer];
            outer_weights.erase(max_outer);
        }
        
        // Move ONLY h_m (S is a fixed model parameter, don't touch it)
        if (max_outer < h_m.size() && empty_index < h_m.size()) {
            h_m[empty_index] = h_m[max_outer];
            h_m[max_outer] = 0;
        }
    }
}

void labels::compact_all() {
    // ========== PART 1: Compact INNER clusters for each outer ==========
    
    for (int mm = 0; mm < M; mm++) {
        // Find all allocated inner for this outer
        std::set<int> allocated_inner;
        for (size_t i = 0; i < inner_allocations.size(); i++) {
            if (outer_allocations[i] == mm) {
                allocated_inner.insert(inner_allocations[i]);
            }
        }
        
        if (allocated_inner.empty()) {
            h_m[mm] = 0;
            continue;
        }
        
        // Create mapping: old_idx → consecutive new_idx (only for allocated clusters)
        std::map<int, int> old_to_new;
        int new_idx = 0;
        for (int old_idx : allocated_inner) {
            old_to_new[old_idx] = new_idx++;
        }
        
        // Update h_m
        h_m[mm] = new_idx;
        
        // Reallocate all observations with new indices
        for (size_t i = 0; i < inner_allocations.size(); i++) {
            if (outer_allocations[i] == mm) {
                inner_allocations[i] = old_to_new[inner_allocations[i]];
            }
        }
        
        // *** CRITICAL: Reorganize inner_sigma maintaining ALL entries 0..S[mm]-1 ***
        if (mm < S.size() && inner_sigma.find(mm) != inner_sigma.end()) {
            auto new_inner_sigma_mm = std::unordered_map<int, std::vector<double>>();
            
            // First: move allocated clusters to new consecutive positions (0, 1, 2, ...)
            for (const auto& [old_idx, new_idx] : old_to_new) {
                if (inner_sigma[mm].find(old_idx) != inner_sigma[mm].end()) {
                    new_inner_sigma_mm[new_idx] = inner_sigma[mm][old_idx];
                }
            }
            
            // Then: fill positions from h_m[mm] to S[mm]-1 with unallocated clusters
            int next_free_idx = h_m[mm];
            for (int old_ss = 0; old_ss < S[mm] && next_free_idx < S[mm]; old_ss++) {
                // If old_ss was not allocated, move it to the next free position
                if (old_to_new.find(old_ss) == old_to_new.end()) {
                    if (inner_sigma[mm].find(old_ss) != inner_sigma[mm].end()) {
                        new_inner_sigma_mm[next_free_idx] = inner_sigma[mm][old_ss];
                    } else {
                        // Create default sigma if it didn't exist
                        int d = inner_sigma[mm].begin()->second.size();
                        new_inner_sigma_mm[next_free_idx] = std::vector<double>(d, 0.5);
                    }
                    next_free_idx++;
                }
            }
            
            // Ensure ALL entries from 0 to S[mm]-1 exist
            for (int ss = 0; ss < S[mm]; ss++) {
                if (new_inner_sigma_mm.find(ss) == new_inner_sigma_mm.end()) {
                    int d = inner_sigma[mm].begin()->second.size();
                    new_inner_sigma_mm[ss] = std::vector<double>(d, 0.5);
                }
            }
            
            inner_sigma[mm] = std::move(new_inner_sigma_mm);
        }
        
        // inner_weights: NO need to reorganize, sample_inner_weights recreates them all
    }
    
    // ========== PART 2: Compact OUTER clusters ==========
    
    std::set<int> allocated_outer;
    for (size_t i = 0; i < outer_allocations.size(); i++) {
        allocated_outer.insert(outer_allocations[i]);
    }
    
    std::map<int, int> outer_old_to_new;
    int new_outer = 0;
    for (int old_outer : allocated_outer) {
        outer_old_to_new[old_outer] = new_outer++;
    }
    
    K = new_outer;
    
    for (size_t i = 0; i < outer_allocations.size(); i++) {
        outer_allocations[i] = outer_old_to_new[outer_allocations[i]];
    }
    
    auto new_outer_center = std::unordered_map<int, std::vector<double>>();
    auto new_outer_weights = std::unordered_map<int, double>();
    auto new_inner_sigma = std::unordered_map<int, std::unordered_map<int, std::vector<double>>>();
    auto new_inner_weights = std::unordered_map<int, std::unordered_map<int, double>>();
    std::vector<int> new_h_m(M, 0);
    
    for (const auto& [old_outer, new_outer_idx] : outer_old_to_new) {
        if (outer_center.find(old_outer) != outer_center.end()) {
            new_outer_center[new_outer_idx] = std::move(outer_center[old_outer]);
        }
        
        if (outer_weights.find(old_outer) != outer_weights.end()) {
            new_outer_weights[new_outer_idx] = outer_weights[old_outer];
        }
        
        if (inner_sigma.find(old_outer) != inner_sigma.end()) {
            new_inner_sigma[new_outer_idx] = std::move(inner_sigma[old_outer]);
        }
        
        if (inner_weights.find(old_outer) != inner_weights.end()) {
            new_inner_weights[new_outer_idx] = std::move(inner_weights[old_outer]);
        }
        
        if (old_outer < h_m.size()) {
            new_h_m[new_outer_idx] = h_m[old_outer];
        }
    }
    
    outer_center = std::move(new_outer_center);
    outer_weights = std::move(new_outer_weights);
    inner_sigma = std::move(new_inner_sigma);
    inner_weights = std::move(new_inner_weights);
    h_m = std::move(new_h_m);
}
