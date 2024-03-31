#pragma once
#include "branching_rules_primitives.hpp"

void get_partial_bit_masks_dfs(int v, const vector<vector<int>> & g_adj, const int & red_vertices_bit_mask, const int & solution_bit_mask, int & partial_bit_mask) {
    partial_bit_mask |=  (1<<v);
    for(int nv : g_adj[v]) {
        //if visited or in the solution set, skip 
        if(partial_bit_mask&(1<<nv) || solution_bit_mask&(1<<nv)) continue;
        get_partial_bit_masks_dfs(nv, g_adj, red_vertices_bit_mask, solution_bit_mask, partial_bit_mask);
    }
}


// return a domination_bit_masks[]
// if domination_bit_masks[i] & (1<<j) == 1, that means we can say solution j dominates solution i 
std::vector<int> get_domination_bit_masks(const _Graph &g, const vector<int> & solutions_bit_masks, const int & red_vertices_bit_mask) {

    //the bist masks of partial set for each solution 
    vector<int> partial_set_bit_masks(solutions_bit_masks.size(), 0);
    vector<vector<int>> g_adj = g.get_adjacent_lists();

    //for every solution compute the partial set
    for(int i = 0; i < solutions_bit_masks.size(); ++i){
        partial_set_bit_masks[i] = 0;
        // for every vertex 
        for(int  j = 0; j < g.n; ++j){
            // such that 
            // 1. not in the red vertex 
            // 2. not already marked in the partial set
            // 3. not in the solution set
            // then dfs the vertices to find all the connected partial vertices
            // if(!(partial_set_bit_masks[i] & (1<<j)) && !(solutions_bit_masks[i] & (1<<j)) && !(red_vertices_bit_mask & (1<<j))){
            if(!((partial_set_bit_masks[i] |  solutions_bit_masks[i] | red_vertices_bit_mask) & (1<<j))){
                get_partial_bit_masks_dfs(j, g_adj, red_vertices_bit_mask, solutions_bit_masks[i], partial_set_bit_masks[i]);
            }
        }
    }


    vector<int> domination_bit_masks(solutions_bit_masks.size(), 0);
    //go over all solution pairs
    for(int i = 0; i < solutions_bit_masks.size(); ++i){
        for(int j = 0; j < i; ++j){
            // if i dominates j 
            if( is_subset_bit_masks(solutions_bit_masks[i], solutions_bit_masks[j]) || 
                size_bit_masks(set_minus_bit_masks(partial_set_bit_masks[i], partial_set_bit_masks[j])) 
            <= size_bit_masks(solutions_bit_masks[j]) - size_bit_masks(solutions_bit_masks[i]) ){
                domination_bit_masks[j] |= (1<<i);
            }
            // if j dominates i
            if( is_subset_bit_masks(solutions_bit_masks[j], solutions_bit_masks[i]) || 
                size_bit_masks(set_minus_bit_masks(partial_set_bit_masks[j], partial_set_bit_masks[i])) 
            <= size_bit_masks(solutions_bit_masks[i]) - size_bit_masks(solutions_bit_masks[j])){
                domination_bit_masks[i] |= (1<<j);
            }
        }
    }

    return domination_bit_masks;
}


pair<vector<vector<int>>,vector<vector<int>>> construct_dominance_graph(const vector<vector<int>> & solutions, const _Graph & g, const vector<int> & red_vertices) {
    int dg_n = solutions.size();
    vector<vector<int>> dg_adj(dg_n);


    // we might define another set to proove, here we return an empty set  
    vector<vector<int>> dg_adj_rv_subset(dg_n);

    // i don't understand why we should use the mapper solution_to_v
    // ??? 
    map<vector<int>, int> solution_to_v;
    for(int i = 0; i < solutions.size(); ++i) {
        solution_to_v[solutions[i]] = i;
    }

    vector<int> solutions_bit_masks = convert_vectors_to_bit_masks(solutions);
    int red_vertices_bit_mask = vector_to_bit_mask(red_vertices);

    vector<int> domination_bit_masks(dg_n);
    domination_bit_masks = get_domination_bit_masks(g, solutions_bit_masks, red_vertices_bit_mask);

    for(int s_i = 0; s_i < solutions.size(); ++s_i) {
        // auto _is_solution_i_dominated_by_solution_j_domination_masks_r = _is_solution_i_dominated_by_solution_j_domination_masks(
            // s_i, solutions_bit_masks, g_adj_lists_bit_masks, safe_red_vertices_subsets_bit_masks);
        // const auto & solution_j_domination_masks = _is_solution_i_dominated_by_solution_j_domination_masks_r.first;
        // const auto & solution_j_rv_subset = _is_solution_i_dominated_by_solution_j_domination_masks_r.second;
        for(int s_j = 0; s_j < solutions.size(); ++ s_j) {
            if(s_i == s_j) continue;
            // if(solution_j_domination_masks[s_j]) {
            if(domination_bit_masks[s_i] & (1<<s_j)){
                dg_adj[solution_to_v[solutions[s_i]]].push_back(solution_to_v[solutions[s_j]]);
                // dg_adj_rv_subset[solution_to_v[solutions[s_i]]].push_back(solution_j_rv_subset[s_j]);
            }
        }
    }

    return {dg_adj, dg_adj_rv_subset};
}

void _strongly_connected_components_dfs1(int v, const vector<vector<int>> & dg_adj, vector<int> & dg_dfs_closed, int & dg_dfs_closed_num) {
    dg_dfs_closed[v] = 0;
    for(int nv : dg_adj[v]) {
        if(dg_dfs_closed[nv]!=-1) continue;
        _strongly_connected_components_dfs1(nv, dg_adj, dg_dfs_closed, dg_dfs_closed_num);
    }
    dg_dfs_closed[v] = dg_dfs_closed_num++;
}

void _strongly_connected_components_dfs2(int v, const vector<vector<int>> & dg_adj, vector<int> & dg_dfs_components, int dg_dfs_components_num) {
    dg_dfs_components[v] = dg_dfs_components_num;
    for(int nv : dg_adj[v]) {
        if(dg_dfs_components[nv]!=-1) continue;
        _strongly_connected_components_dfs2(nv, dg_adj, dg_dfs_components, dg_dfs_components_num);
    }
}

std::pair<vector<vector<int>>, pair<vector<vector<int>>, vector<vector<int>>>> compute_all_dominance_free_solutions(const vector<vector<int>> & solutions, const _Graph & g, const vector<int> & red_vertices) {
    int dg_n = solutions.size();
    map<vector<int>, int> solution_to_v;
    for(int i = 0; i < solutions.size(); ++i) {
        solution_to_v[solutions[i]] = i;
    }

    //construc the dominance graph 
    //here the second element of the pair is useless
    auto construct_dominance_graph_r = construct_dominance_graph(solutions, g, red_vertices);
    const auto & dg_adj = construct_dominance_graph_r.first;
    const auto & dg_adj_rv_subset = construct_dominance_graph_r.second;

    //dg_dfs_close[i] means the order of i in the dfs sequence 
    vector<int> dg_dfs_closed(dg_n, -1);
    int dg_dfs_closed_num = 0;
    for(int i = 0; i < dg_n; ++i) {
        if(dg_dfs_closed[i]==-1) {
            _strongly_connected_components_dfs1(i, dg_adj, dg_dfs_closed, dg_dfs_closed_num);
        }
    }

    //the reversed dominance graph 
    vector<vector<int>> dg_adj_reversed(dg_n);
    for(int i = 0; i < dg_n; ++i) {
        for(int nv : dg_adj[i]) {
            dg_adj_reversed[nv].push_back(i);
        }
    }

    //if there is an edge (i,j) in the dominance graph
    // that means solution i is dominated by solution j 
    // thus we have to dfs2 in the reversed order of dg_dfs_close
    vector<int> dg_dfs_closed_order(dg_n);
    for(int i = 0; i < dg_n; ++i) {
        dg_dfs_closed_order[dg_n-1-dg_dfs_closed[i]] = i;
    }
    // the num i vertex in dominance graph belongs to the component dg_dfs_components[i]
    vector<int> dg_dfs_components(dg_n, -1);
    // the number of strongly components in dominance graph
    int dg_dfs_components_num = 0;
    for(int v : dg_dfs_closed_order) {
        if(dg_dfs_components[v]==-1) {
            _strongly_connected_components_dfs2(v, dg_adj_reversed, dg_dfs_components, dg_dfs_components_num++);
        }
    }

    //build the contracted dominance graph  
    int cdg_n = dg_dfs_components_num;
    vector<vector<int>> cdg_adj(cdg_n);
    std::set<pair<int,int>> contraction_edges;
    for(int i = 0; i < dg_n; ++i) {
        for(int nv : dg_adj[i]) {
            int ci = dg_dfs_components[i];
            int cj = dg_dfs_components[nv];
            if(ci==cj)continue; // if i and nv are in the same strongly connected component 
            if(!contraction_edges.count({ci, cj})) {
                contraction_edges.insert({ci, cj});
                cdg_adj[ci].push_back(cj);
                // if there is and edge (ci,cj)
                // that means the solutions of component ci is dominated by the solutions of compnent cj  
            }
        }
    }

    vector<vector<int>> csink_solutions;

    //go over every vertex of contracted dominance graph 
    for(int i = 0; i < cdg_n; ++i) {
        // if i is not dominated by any solutions in other component
        if(cdg_adj[i].size()==0) {
            //choose any solution in this strongly connected component
            for(int j = 0; j < dg_n; ++j) {
                if(dg_dfs_components[j]==i) {
                    csink_solutions.push_back(solutions[j]);
                    break; // notice here break, any one is OK
                }
            }
        }
    }

    return {csink_solutions, {dg_adj, dg_adj_rv_subset}};
}


