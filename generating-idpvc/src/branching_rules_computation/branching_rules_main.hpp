#pragma once
#include "branching_rules_primitives.hpp"
#include "reduction_rules.hpp"
#include "minimal_solutions_computation.hpp"
#include "dominance_free_solutions_computation.hpp"
#include "adjusted_solutions_computation.hpp"

BranchingRule generate_branching_rule(const _Graph & g) {
    vector<int> red_vertices;
    for(int i = 0; i < g.n; ++i) {
        if(g.red_vertices_mask & (1<<i)) {
            red_vertices.push_back(i);
        }
    }

    #ifndef NOREDUCTION_RULES

    if(is_g_rv_handled_by_red_component_reduction(g, red_vertices)) {
        return BranchingRule(g, red_vertices, {}, {}, {}, -1, "handled_rcr");
    }

    if(is_g_rv_handled_by_red_star_reduction(g, red_vertices)) {
        return BranchingRule(g, red_vertices, {}, {}, {}, -1, "handled_rsr");
    }

    if(is_g_rv_handled_by_vertex_cover_struction(g, red_vertices)) {
        return BranchingRule(g, red_vertices, {}, {}, {}, -1, "handled_vc_struction");
    }

    if(is_g_rv_handled_by_vertex_cover_dominance(g, red_vertices)) {
        return BranchingRule(g, red_vertices, {}, {}, {}, -1, "handled_vc_dominance");
    }

    #endif

    const auto & minimal_solutions = compute_all_solutions(g);
    // const auto & minimal_solutions = compute_all_minimal_solutions(g);
    double bf_minimal_solutions = compute_branching_factor_for_solutions(minimal_solutions);

    if(bf_minimal_solutions <= DPVC_BF)
    {
        proof_branching_rule(g, minimal_solutions, {}, {}, {});
        return BranchingRule(g, {}, minimal_solutions, {}, minimal_solutions, bf_minimal_solutions, "minimal");
    }
    #ifdef NODOMINANCE
    return BranchingRule(g, {}, minimal_solutions, {}, minimal_solutions, bf_minimal_solutions, "minimal");
    #endif

    const auto & compute_all_dominance_free_solutions_r = compute_all_dominance_free_solutions(minimal_solutions, g, red_vertices);
    const auto & dominance_free_solutions = compute_all_dominance_free_solutions_r.first;
    double bf_dominance_free_solutions = compute_branching_factor_for_solutions(dominance_free_solutions);

    if(bf_dominance_free_solutions <= DPVC_BF)
    {
        proof_branching_rule(g, minimal_solutions, dominance_free_solutions, compute_all_dominance_free_solutions_r.second, dominance_free_solutions);
        return BranchingRule(g, red_vertices, minimal_solutions, dominance_free_solutions, dominance_free_solutions, bf_dominance_free_solutions, "dominance_free");
    }
    #ifdef NOADJUSTMENT
    return BranchingRule(g, red_vertices, minimal_solutions, dominance_free_solutions, dominance_free_solutions, bf_dominance_free_solutions, "dominance_free");
    #endif

    const auto & adjusted_solutions = compute_adjusted_solutions_bit_masks(g, dominance_free_solutions);
    double bf_adjusted_solutions = compute_branching_factor_for_solutions(adjusted_solutions);

    if(bf_adjusted_solutions <= DPVC_BF) {
        proof_branching_rule(g, minimal_solutions, dominance_free_solutions, compute_all_dominance_free_solutions_r.second, adjusted_solutions);
    }

    return BranchingRule(g, red_vertices, minimal_solutions, dominance_free_solutions, adjusted_solutions, bf_adjusted_solutions, "adjusted_dominance");
}
