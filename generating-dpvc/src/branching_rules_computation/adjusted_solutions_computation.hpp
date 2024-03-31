#pragma once
#include "branching_rules_primitives.hpp"


vector<vector<int>> compute_adjusted_solutions(const _Graph & g, const vector<vector<int>> & solutions) {
    assert(DPVC_BF>0);

    vector<vector<int>> current_solutions = solutions;
    double bf_current_solutions = compute_branching_factor_for_solutions(current_solutions);

    std::set<std::vector<int>> adjusting_branches_set;
    for(const auto & branch : current_solutions) {
        for(const auto & adjusting_branch : ordered_powerset_nonempty(branch)) adjusting_branches_set.insert(adjusting_branch);
    }

    std::vector<std::vector<int>> neg_diff_adjusting_branches;
    for(const auto & adjusting_branch : adjusting_branches_set) {
        double difference = pow(DPVC_BF, -(int)adjusting_branch.size());
        for(const auto & branch : current_solutions) {
            if(is_subset(adjusting_branch, branch)) difference -= pow(DPVC_BF, -(int)branch.size());
        }
        if(difference<0)neg_diff_adjusting_branches.push_back(adjusting_branch);
    }

    bool something_changed = true;

    while(something_changed) {
        something_changed = false;

        vector<int> best_av_subset;
        double best_av_bf=-1;

        for(const auto & neg_diff_adjusting_branch : neg_diff_adjusting_branches) {
            vector<int> adjusted_branching_vector;
            adjusted_branching_vector.push_back(neg_diff_adjusting_branch.size());
            for(const vector<int> & solution : current_solutions) {
                if(is_subset(neg_diff_adjusting_branch, solution) == false) {
                    adjusted_branching_vector.push_back(solution.size());
                }
            }
            double adjusted_branching_factor = compute_branching_factor_for_branching_vector(adjusted_branching_vector);
            if(best_av_bf == -1 || adjusted_branching_factor < best_av_bf) {
                best_av_bf = adjusted_branching_factor;
                best_av_subset = neg_diff_adjusting_branch;
            }
        }

        if(best_av_bf != -1 && best_av_bf < bf_current_solutions) {
            vector<vector<int>> new_current_solutions;
            new_current_solutions.push_back(best_av_subset);
            for(const vector<int> & solution : current_solutions) {
                if(is_subset(best_av_subset, solution) == false) {
                    new_current_solutions.push_back(solution);
                }
            }
            swap(current_solutions, new_current_solutions);
            bf_current_solutions = best_av_bf;
            something_changed = true;

        }
    }

    return current_solutions;
}


vector<vector<int>> compute_adjusted_solutions_bit_masks(const _Graph & g, const vector<vector<int>> & solutions) {
    assert(DPVC_BF>0);

    vector<vector<int>> current_solutions = solutions;
    double bf_current_solutions = compute_branching_factor_for_solutions(current_solutions);

    // We construct the possible adjustments we can do.
    std::set<int> adjusting_branches_set;
    for(const auto & branch : current_solutions) {
        for(const auto & adjusting_branch : ordered_powerset_nonempty(branch)) adjusting_branches_set.insert(vector_to_bit_mask(adjusting_branch));
    }

    // OPTIMIZATION:
    // We can filter out adjustments that have no chance to improve the branching factor.
    vector<int> neg_diff_adjusting_branches_bit_masks;
    vector<int> current_solutions_bit_masks = convert_vectors_to_bit_masks(current_solutions);

    for(int adjusting_branch : adjusting_branches_set) {
        double difference = pow(DPVC_BF, -(int)size_bit_masks(adjusting_branch));
        for(int branch : current_solutions_bit_masks) {
            if(is_subset_bit_masks(adjusting_branch, branch)) difference -= pow(DPVC_BF, -(int)size_bit_masks(branch));
        }
        if(difference<0)neg_diff_adjusting_branches_bit_masks.push_back(adjusting_branch);
    }

    // Sometimes trying to cover everything by two branches works.
    for(int i = 0; i < neg_diff_adjusting_branches_bit_masks.size(); ++i) {
        for(int j = i+1; j < neg_diff_adjusting_branches_bit_masks.size(); ++j) {
            int branch1 = neg_diff_adjusting_branches_bit_masks[i];
            int branch2 = neg_diff_adjusting_branches_bit_masks[j];
            double bf = compute_branching_factor_for_branching_vector({size_bit_masks(branch1), size_bit_masks(branch2)});
            if(bf>DPVC_BF) continue;

            bool covers = true;
            for(int solution : current_solutions_bit_masks) {
                if(is_subset_bit_masks(branch1, solution) == false && is_subset_bit_masks(branch2, solution) == false) {
                    covers = false;
                    break;
                }
            }
            if(covers) {
                return convert_bit_masks_to_vectors({branch1, branch2});
            }
        }
    }

    // Repeat until nothing changes - no improvement can be made.
    bool something_changed = true;
    while(something_changed) {
        something_changed = false;

        int best_av_subset;
        double best_av_bf=-1;

        for(int neg_diff_adjusting_branch : neg_diff_adjusting_branches_bit_masks) {
            vector<int> adjusted_branching_vector;
            adjusted_branching_vector.push_back(size_bit_masks(neg_diff_adjusting_branch));
            for(int solution : current_solutions_bit_masks) {
                if(is_subset_bit_masks(neg_diff_adjusting_branch, solution) == false) {
                    adjusted_branching_vector.push_back(size_bit_masks(solution));
                }
            }
            double adjusted_branching_factor = compute_branching_factor_for_branching_vector(adjusted_branching_vector);
            if(best_av_bf == -1 || adjusted_branching_factor < best_av_bf) {
                best_av_bf = adjusted_branching_factor;
                best_av_subset = neg_diff_adjusting_branch;
            }
        }

        if(best_av_bf != -1 && best_av_bf < bf_current_solutions) {
            vector<int> new_current_solutions_bit_masks;
            new_current_solutions_bit_masks.push_back(best_av_subset);
            for(int solution : current_solutions_bit_masks) {
                if(is_subset_bit_masks(best_av_subset, solution) == false) {
                    new_current_solutions_bit_masks.push_back(solution);
                }
            }
            swap(current_solutions_bit_masks, new_current_solutions_bit_masks);
            bf_current_solutions = best_av_bf;
            something_changed = true;

        }
    }

    return convert_bit_masks_to_vectors(current_solutions_bit_masks);
}