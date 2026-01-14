#pragma once

#include <unordered_map>
#include <vector>
#include <set>
#include <algorithm>

class labels {
private:

    std::vector<int> inner_allocations;
    std::vector<int> outer_allocations;

    std::unordered_map<int, std::unordered_map<int, std::vector<double>>> inner_sigma; // chiave esterna: outer cluster, 
    //chiave interna: inner cluster locale, valore: vettore di sigma per ogni dimensione
    std::unordered_map<int, std::vector<double>> outer_center; // chiave: outer cluster, valore: vettore di centri per ogni variabile

    std::unordered_map<int, double> outer_weights; // chiave: outer cluster, valore: peso w_un
    std::unordered_map<int, std::unordered_map<int, double>> inner_weights; // chiave esterna: outer cluster, chiave interna: inner cluster locale, valore: peso q_un

    int M = 0; // Numero di cluster outer totali
    int K = 0; // Numero di cluster outer allocati
    std::vector<int> h_m; // Numero di cluster inner allocati per ogni outer cluster 
    std::vector<int> S; // Numero di cluster inner totali per ogni outer cluster

    void compact_allocations(const int empty_index, bool inner = true, const int outer_index = -1); // Funzione privata per compattare le allocazioni

public:
    // Default constructor
    labels() : M(0), K(0) {}

    // Parameterized constructor
    labels(const std::vector<int>& inner_allocs,
           const std::vector<int>& outer_allocs,
           const std::unordered_map<int, std::unordered_map<int, std::vector<double>>>& inner_sigma_,
           const std::unordered_map<int, std::vector<double>>& outer_center_,
           const std::unordered_map<int, double>& outer_weights_,
           const std::unordered_map<int, std::unordered_map<int, double>>& inner_weights_,
           int M_,
           const std::vector<int>& S_)
        : inner_allocations(inner_allocs),
          outer_allocations(outer_allocs),
          inner_sigma(inner_sigma_),
          outer_center(outer_center_),
          outer_weights(outer_weights_),
          inner_weights(inner_weights_),
          M(M_),
          S(S_) {
        
        std::set<int> unique_m(outer_allocations.begin(), outer_allocations.end());
        K = unique_m.size();

        h_m.resize(M, 0);
        for (int i = 0; i < inner_allocations.size(); i++) {
            int mm = outer_allocations[i];
            int ss = inner_allocations[i];
            // Aggiorna il massimo visto per questo outer cluster
            if (ss + 1 > h_m[mm]) {
                h_m[mm] = ss + 1;
            }
        }   
    }

    // ========== Getters ==========
    const std::vector<int>& get_inner_allocations() const {
        return inner_allocations;
    }

    const std::vector<int>& get_outer_allocations() const {
        return outer_allocations;
    }

    int get_inner_allocation_at(int index) const {
        return inner_allocations.at(index);
    }

    int get_outer_allocation_at(int index) const {
        return outer_allocations.at(index);
    }

    const std::unordered_map<int, std::unordered_map<int, std::vector<double>>>& get_inner_sigma() const {
        return inner_sigma;
    }

    const std::unordered_map<int, std::vector<double>>& get_outer_center() const {
        return outer_center;
    }

    const std::unordered_map<int, double>& get_outer_weights() const {
        return outer_weights;
    }

    const std::unordered_map<int, std::unordered_map<int, double>>& get_inner_weights() const {
        return inner_weights;
    }

    int get_M() const {
        return M;
    }

    int get_K() const {
        return K;
    }

    const std::vector<int>& get_h_m() const {
        return h_m;
    }

    const std::vector<int>& get_S() const {
        return S;
    }


    // ========== Setters ==========
    void set_inner_allocations(std::vector<int>& inner_allocs) {
        inner_allocations = std::move(inner_allocs);
    }

    void set_outer_allocations(std::vector<int>& outer_allocs) {
        outer_allocations = std::move(outer_allocs);
    }

    void set_inner_allocation_at(int index, int value) {
        int old_inner_idx = inner_allocations.at(index);
        int outer_idx = outer_allocations.at(index);
        inner_allocations.at(index) = value;

        // Check if old_inner_idx is still occupied
        bool still_occupied = false;
        for (size_t i = 0; i < inner_allocations.size(); i++) {
            if (outer_allocations[i] == outer_idx && inner_allocations[i] == old_inner_idx) {
                still_occupied = true;
                break;
            }
        }

        if (!still_occupied) {
            compact_allocations(old_inner_idx, true, outer_idx);
        }
    }

    void set_outer_allocation_at(int index, int value) {
        int old_outer_idx = outer_allocations.at(index);
        outer_allocations.at(index) = value;

        // Check if old_outer_idx is still occupied
        auto it = std::find(outer_allocations.begin(), outer_allocations.end(), old_outer_idx);
        if(it == outer_allocations.end()) {
            compact_allocations(old_outer_idx, false);
        }
    }

    void set_inner_sigma(const std::unordered_map<int, std::unordered_map<int, std::vector<double>>>& inner_sigma_) {
        inner_sigma = inner_sigma_;
    }

    void set_outer_center(const std::unordered_map<int, std::vector<double>>& outer_center_) {
        outer_center = outer_center_;
    }

    void set_outer_weights(const std::unordered_map<int, double>& outer_weights_) {
        outer_weights = outer_weights_;
    }

    void set_inner_weights(const std::unordered_map<int, std::unordered_map<int, double>>& inner_weights_) {
        inner_weights = inner_weights_;
    }

    void set_M(int M_) {
        M = M_;
    }

    void set_h_m(const std::vector<int>& h_m_) {
        h_m = h_m_;
    }

    void set_S(const std::vector<int>& S_) {
        S = S_;
    }
};