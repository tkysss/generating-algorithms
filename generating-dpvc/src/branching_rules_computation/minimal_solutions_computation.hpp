#pragma once
#include "branching_rules_primitives.hpp"

void _compute_all_minimal_solutions_dfs(const _Graph & g, int cur_v, vector<int> & cur_sol, vector<vector<int>> & out_solutions) {
    cur_sol.push_back(cur_v);

    if(is_g_solved_by_solution(g, cur_sol)) {
        out_solutions.push_back(cur_sol);
    }
    else {
        for(int next_v = cur_v + 1; next_v < g.n; ++next_v) {
            _compute_all_minimal_solutions_dfs(g, next_v, cur_sol, out_solutions);
        }
    }
    cur_sol.pop_back();
}

vector<vector<int>> compute_all_minimal_solutions(const _Graph & g) {
    vector<vector<int>> solutions;
    vector<int> cur_sol;
    for(int i = 0; i < g.n; ++i) {
        _compute_all_minimal_solutions_dfs(g, i, cur_sol, solutions);
    }

    vector<vector<int>> f_solutions;

    for(int idx_i = 0; idx_i < solutions.size(); ++idx_i) {
        bool ok = true;
        for(int idx_j = idx_i + 1; idx_j < solutions.size(); ++idx_j) {
            if(is_subset(solutions[idx_j], solutions[idx_i])) {
                ok = false;
                break;
            }
        }
        if(ok){
            f_solutions.push_back(solutions[idx_i]);
        }
    }

    return f_solutions;
}
