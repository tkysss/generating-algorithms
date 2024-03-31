#pragma once
#include "../graph_computation/graph_primitives.hpp"

struct BranchingRule {
    _Graph g;
    vector<int> red_vertices;
    vector<vector<int>> minimal_branches;
    vector<vector<int>> dominance_free_branches;
    vector<vector<int>> adjusted_branches;
    double bf;
    std::string type;

    BranchingRule() : g(0,_GraphEdges(), make_empty_nauty_canonical_rep_with_relabeling(),0) {}

    BranchingRule(const _Graph & g,
    const vector<int> & red_vertices,
    const vector<vector<int>> & minimal_branches,
    const vector<vector<int>> & dominance_free_branches,
    const vector<vector<int>> & adjusted_branches,
    double bf,
    std::string type
    ) :
    g(g),
    red_vertices(red_vertices),
    minimal_branches(minimal_branches),
    dominance_free_branches(dominance_free_branches),
    adjusted_branches(adjusted_branches),
    bf(bf),
    type(type)
    {}

};

bool is_g_solved_by_solution(const _Graph & g, const vector<int> & solution) {
    assert(DPVC_PATH_LEN>0);
    if(g.n - solution.size() < DPVC_PATH_LEN) return true;

    if(DPVC_PATH_LEN == 2) {
        vector<int> v_map(g.n, 0);
        for(int v : solution)v_map[v] = 1;
        for(auto e : g.get_edges_list()) {
            if(!(v_map[e.first] || v_map[e.second]))return false;
        }
        return true;
    }

    // if(DPVC_PATH_LEN == 3) {
    //     vector<int> v_map(g.n, 0);
    //     for(int v : solution)v_map[v] = 1;
    //     vector<int> v_degs(g.n, 0);
    //     for(auto e : g.get_edges_list()) {
    //         if(v_map[e.first] || v_map[e.second])continue;
    //         v_degs[e.first]++;
    //         v_degs[e.second]++;
    //         if(v_degs[e.first] >= 2 || v_degs[e.second] >= 2)return false;
    //     }

    //     return true;
    // }

    for(vector<int> path_vertices : ordered_combinations(set_minus(range(g.n), solution), DPVC_PATH_LEN)) {
        int edge_count_on_path_vertices = 0;
        for(int i = 0; i < DPVC_PATH_LEN; ++i){
            for(int j = i+1; j < DPVC_PATH_LEN; ++j){
                int a = path_vertices[i];
                int b = path_vertices[j];
                if(a>b)std::swap(a,b);
                edge_count_on_path_vertices += g.has_edge(a,b);
            }
        }
        // Surely no path on these if there are not enough edges
        // if(edge_count_on_path_vertices < DPVC_PATH_LEN - 1) continue;
        if(edge_count_on_path_vertices != DPVC_PATH_LEN - 1) continue;

        // Brute force the path.
        do {
            bool has_path = 1;
            for(int i = 0; i < DPVC_PATH_LEN-1; ++i) {
                int a = path_vertices[i];
                int b = path_vertices[i+1];
                if(a>b)std::swap(a,b);
                if(!g.has_edge(a,b)){has_path=0;break;}
            }
            if(has_path) {
                return false;
            }
        } while(next_permutation(path_vertices.begin(), path_vertices.end()));
    }
    return true;
}


int int_pow(int base, int exponent) {
    int r = 1;
    for(int i = 0; i < exponent; ++i) r*=base;
    return r;
}

double compute_branching_factor_for_branching_vector(const vector<int> & branching_vector) {
    // try some int solutions first
    if(DPVC_BF >=2) {
        int mx = 0;
        for(int b : branching_vector)if(b>mx)mx=b;
        for(int mi=2;mi<=DPVC_BF; ++mi) {
            int v = -int_pow(mi,  mx);
            for(int b : branching_vector) v+= int_pow(mi, mx-b);
            if(v==0) return mi;
        }
    }

    double lo = 1;
    double hi = DPVC_PATH_LEN+1;
    const double ACCURACY = 100000;
    const double INV_ACCURACY = 1./ACCURACY;
    while(hi - lo > INV_ACCURACY) {
        double mi = (hi + lo) / 2.;
        double v = -1;
        for(int b : branching_vector)
            v += pow(mi, -1.*b);
        if(v > 0)lo = mi;
        else hi = mi;
    }
    return ceil(hi * ACCURACY) * INV_ACCURACY;
}

double compute_branching_factor_for_solutions(const vector<vector<int>> & solutions) {
    vector<int> branching_vector;
    for(const vector<int> & solution : solutions)branching_vector.push_back(solution.size());
    return compute_branching_factor_for_branching_vector(branching_vector);
}
