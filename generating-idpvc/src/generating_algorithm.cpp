#include <chrono>
#include "common_primitives/primitive_operations.hpp"
#include "common_primitives/output_operations.hpp"
#include "graph_computation/graph_operations.hpp"
#include "branching_rules_computation/branching_rules_main.hpp"

void _set_global_variables(int argc, char ** argv) {
    if(argc == 4) {
        DPVC_PATH_LEN = std::stoi(argv[1]);
        DPVC_BF = std::stod(argv[2]);
        GENERATIONS = std::stoi(argv[3]);
    }
    else {
        cerr<<"run like ./generating_algorithm DPVC_PATH_LEN DPVC_BF GENERATIONS"<<endl;
        cerr<<"  DPVC_PATH_LEN - integer - the length of d-path"<<endl;
        cerr<<"  DPVC_BF       - float   - the desired branching factor"<<endl;
        cerr<<"  GENERATIONS   - integer - number of generations to run before giving up"<<endl;
        exit(1);
    }
    NUM_THREADS = omp_get_max_threads();
    omp_set_num_threads(NUM_THREADS);
    omp_set_nested(1);

    OUTPUT_FILES_PREFIX = std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]);

    cerr<<"DPVC_PATH_LEN: "<<DPVC_PATH_LEN<<endl;
    cerr<<"DPVC_BF: "<<DPVC_BF<<endl;
    cerr<<"GENERATIONS: "<<GENERATIONS<<endl;
    cerr<<"NUM_THREADS: "<<NUM_THREADS<<endl;
}


// Returns the set of initial bumpy graphs - set of all connected non-isomorphic bumpy graphs with DPVC_PATH_LEN vertices.
vector<_Graph> get_initial_bumpy_graphs(vector<_Nauty> & nts) {
    vector<_Graph> expanded_graphs;
    _ForbiddenInducedSubgraphs fis;
    expanded_graphs.push_back(_Graph(1, {}, {}, nts[0]));
    for(int i = 0; i < DPVC_PATH_LEN - 1; ++i) {
        expanded_graphs = gather_nonisomorphic_graphs(parallel_expand_graph_without_forbidden_subgraph(expanded_graphs, fis, nts));
    }
    vector<_Graph> current_bumpy_graphs;
    for(const _Graph & g : expanded_graphs) {
        if(!is_g_solved_by_solution(g, {})) current_bumpy_graphs.push_back(g);
    }

    return current_bumpy_graphs;
}

// Returns the set of initial forbidden induced subgraphs. In the start, it is empty.
_ForbiddenInducedSubgraphs get_initial_forbidden_induced_subgraphs(vector<_Nauty> & nts) {
    return _ForbiddenInducedSubgraphs();
}

int main(int argc, char ** argv) {
    _set_global_variables(argc, argv);
    proof_open_proof_file();
    proof_parameters(DPVC_PATH_LEN, DPVC_BF, GENERATIONS);

    std::chrono::high_resolution_clock::time_point _time_start = std::chrono::high_resolution_clock::now();

    vector<_Nauty> nts(NUM_THREADS);

    vector<BranchingRule> branching_rules;
    vector<BranchingRule> not_solved_branching_rules;

    _ForbiddenInducedSubgraphs forbidden_induced_subgraphs = get_initial_forbidden_induced_subgraphs(nts);
    vector<_Graph> current_bumpy_graphs = get_initial_bumpy_graphs(nts);

    // We start from generation 1 and repeat the process for GENERATIONS unless the algorithm is found.
    for(int generation = 1; generation <= GENERATIONS; ++generation) {
        cerr<<"-- generation="<<generation<<endl;
        cerr<<"-- starting with current_bumpy_graphs="<<current_bumpy_graphs.size()<<endl;

        proof_generation_start(generation, current_bumpy_graphs);

        // The generation consists of rounds.
        // We start from round 1 and we repeat the rounds until nothing changes (no new branching rules with good enough branching factor).
        int round = 1;
        bool something_changed = true;

        // In the first round, we expand all current_bumpy_graphs agains the global forbidden_induced_subgraphs.
        proof_expansion_start();
        // OPTIMIZATION:
        // After that we need to filter the expansions only against the local_forbidden_induced_subgraphs as the global ones were already tested.
        vector<vector<_Graph>> current_bumpy_graph_expansions = parallel_expand_graph_without_forbidden_subgraph(
            current_bumpy_graphs, forbidden_induced_subgraphs, nts);

        proof_expansion_end();

        _ForbiddenInducedSubgraphs local_forbidden_induced_subgraphs;

        // For each graph in current_bumpy_graphs, we will generate a branching rule.
        vector<BranchingRule> current_brs(current_bumpy_graphs.size());

        // OPTIMIZATION:
        // We need to regenerate branching rule for a graph only if its red vertices changed. Otherwise, no new information was gained.
        // This keeps track of the changes.
        vector<uint32_t> red_vertices_changed_masks(current_bumpy_graphs.size());

        while(something_changed) {
            cerr<<"-- round="<<round<<" with current_bumpy_graphs="<<current_bumpy_graphs.size()<<endl;

            proof_round_start(round);

            // We make an assumption that nothing will change in this round.
            something_changed = false;

            // If this is not the first round, filter the graphs againts local_forbidden_induced_subgraphs
            cerr<<"-- round="<<round<<" before parallel filter"<<endl;
            if(round!=1) {
                parallel_filter_graph_expansions_containing_forbidden_induced_subgraph(
                    current_bumpy_graphs, current_bumpy_graph_expansions, local_forbidden_induced_subgraphs, nts);
            }
            cerr<<"-- round="<<round<<" after parallel filter"<<endl;

            // Compute the new annotations of red_vertices and keep track for which graphs the red vertices changed.
            #pragma omp parallel default(none) shared(current_bumpy_graphs, current_bumpy_graph_expansions, red_vertices_changed_masks)
            {
                #pragma omp for schedule(dynamic)
                for(int g_idx=0; g_idx < current_bumpy_graphs.size(); ++g_idx) {
                    red_vertices_changed_masks[g_idx] = current_bumpy_graphs[g_idx].annotate_red_vertices(current_bumpy_graph_expansions[g_idx]);
                }
            }

            // Gather the indexes of the graphs with changed red vertices
            // In the first round, the annotation must be recomputed for all, so we consider all of them as changed.
            vector<int> g_idx_with_changed_red_vertices_masks;
            for(int g_idx=0; g_idx < current_bumpy_graphs.size(); ++g_idx) {
                if(round==1 || red_vertices_changed_masks[g_idx]) g_idx_with_changed_red_vertices_masks.push_back(g_idx);
            }

            cerr<<"-- round="<<round<<" after annotation "<<g_idx_with_changed_red_vertices_masks.size()<<endl;


            int _progress = 0;
            int _progress_s = (g_idx_with_changed_red_vertices_masks.size()+100-1) / 100;
            auto _progress_clock = std::chrono::steady_clock::now();
            int64_t _progress_duration = 0;
            cerr<<"-- generating rules------------------------------------->"<<g_idx_with_changed_red_vertices_masks.size()<<endl;

            // We generate new branching rules for only those graphs whose red vertices changed.
            #pragma omp parallel default(none) shared(current_bumpy_graphs, current_brs, g_idx_with_changed_red_vertices_masks, cerr, _progress, _progress_s, _progress_clock, _progress_duration)
            {
                #pragma omp for schedule(dynamic)
                for(int gg_idx=0; gg_idx < g_idx_with_changed_red_vertices_masks.size(); ++gg_idx) {
                    int g_idx = g_idx_with_changed_red_vertices_masks[gg_idx];
                    current_brs[g_idx] = generate_branching_rule(current_bumpy_graphs[g_idx]);

                    #pragma omp atomic update
                    _progress ++;
                    if((_progress % _progress_s) == 0){
                        #pragma omp critical
                        {
                        _progress_duration += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - _progress_clock).count();
                        int64_t eta_nanos = (g_idx_with_changed_red_vertices_masks.size() - _progress) * (_progress_duration / _progress);
                        cerr<<eta_nanos/1000/1000/1000<<"s ";
                        _progress_clock = std::chrono::steady_clock::now();
                        }
                    }
                }
            }
            cerr<<endl;

            cerr<<"-- round="<<round<<" after generating"<<endl;

            // We go over all generated branching rules and decide if they are good or bad.
            for(int g_idx=0; g_idx < current_bumpy_graphs.size(); ++g_idx) {
                const _Graph & g = current_bumpy_graphs[g_idx];
                const BranchingRule & br = current_brs[g_idx];
                // If the branching rule is good, something changed and we move it from current branching rules to good branching rules.
                // And we adjust current_bumpy_graphs, forbidden_induced_subgraphs, ... accordingly.
                if(br.bf <= DPVC_BF) {
                    something_changed = true;

                    forbidden_induced_subgraphs.add(g.canonical_rep_with_relabeling.canonical_rep, g.n);
                    local_forbidden_induced_subgraphs.add(g.canonical_rep_with_relabeling.canonical_rep, g.n);
                    branching_rules.push_back(br);

                    // proof_new_branching_rule(g);
 
                    current_bumpy_graphs[g_idx] = current_bumpy_graphs[current_bumpy_graphs.size()-1];
                    current_bumpy_graphs.pop_back();
                    current_bumpy_graph_expansions[g_idx] = current_bumpy_graph_expansions[current_bumpy_graph_expansions.size()-1];
                    current_bumpy_graph_expansions.pop_back();
                    current_brs[g_idx] = current_brs[current_brs.size()-1];
                    current_brs.pop_back();
                    assert(current_bumpy_graphs.size() == current_bumpy_graph_expansions.size());
                    assert(current_bumpy_graphs.size() == current_brs.size());
                    --g_idx;
                }
                else {
                    // Otherwise the branching rule is bad and we leave it be for another round.
                }
            }

            cerr<<"-- round="<<round<<" after decision"<<endl;

            proof_round_end(round);

            // If we happen to solve all current_bumpy_graphs, the algorithm is complete and we end!
            if(current_bumpy_graphs.size() == 0) {
                cerr<<"-- current_bumpy_graphs empty"<<endl;
                break;
            }
            // Otherwise we go to the next round.
            round++;

            // If this is the last generation, store the computed branching_rules, that werent solved.
            if(something_changed == false && generation == GENERATIONS) {
                not_solved_branching_rules = current_brs;
            }
        }

        // OPTIMIZATION:
        // For each graph, propagate its red_vertices into its expansions. As the expansion is a supergraph, it must have the at least the same set of red_vertices.
        for(int g_idx=0; g_idx < current_bumpy_graphs.size(); ++g_idx) {
            for(_Graph & eg : current_bumpy_graph_expansions[g_idx]) {
                eg.red_vertices_mask = current_bumpy_graphs[g_idx].red_vertices_mask;
            }
        }

        proof_generation_end(generation, current_bumpy_graphs, current_bumpy_graph_expansions);

        // The generation ended. If this was the last iteration, dont bother with gathering.
        if(generation==GENERATIONS || current_bumpy_graphs.size() == 0) {
            break;
        }

        // Finally, we gather all non-isomorphic expansions and treat them as new current_bumpy_graphs for the next generation.
        current_bumpy_graphs = gather_nonisomorphic_graphs(current_bumpy_graph_expansions);

        proof_generation_next(generation, current_bumpy_graphs);
    }


    // Final output - find maximum branch factor, output all the rules into their respective files, gather some statistics, and end.
    std::chrono::high_resolution_clock::time_point _time_end = std::chrono::high_resolution_clock::now();
    int _time_s = std::chrono::duration_cast<std::chrono::seconds>(_time_end - _time_start).count();

    double max_branching_rules_bf = 0;
    int max_branching_rules_size = 0;
    for(const BranchingRule & br : branching_rules){
        max_branching_rules_bf = std::max(max_branching_rules_bf, br.bf);
        max_branching_rules_size = std::max(max_branching_rules_size, br.g.n);
    }
    double max_not_solved_branching_rules_bf = 0;
    for(const BranchingRule & br : not_solved_branching_rules){max_not_solved_branching_rules_bf = std::max(max_not_solved_branching_rules_bf, br.bf);}

    if(not_solved_branching_rules.size()==0) {cerr<<"success with:"<<endl;}
    else {cerr<<"failed with:"<<endl;}
    cerr<<"-- max_branching_rules_bf="<<max_branching_rules_bf<<endl;
    cerr<<"-- branching_rules="<<branching_rules.size()<<endl;
    cerr<<"-- max_not_solved_branching_rules_bf="<<max_not_solved_branching_rules_bf<<endl;
    cerr<<"-- not_solved_branching_rules="<<not_solved_branching_rules.size()<<endl;

    std::ofstream stat_br_ofs (OUTPUT_FILES_PREFIX + ".stats.out", std::ofstream::out);
    stat_br_ofs<<"{"<<endl;
    stat_br_ofs<<"\"max_branching_rules_bf\":"<<max_branching_rules_bf<<","<<endl;
    stat_br_ofs<<"\"branching_rules\":"<<branching_rules.size()<<","<<endl;
    stat_br_ofs<<"\"max_branching_rules_size\":"<<max_branching_rules_size<<","<<endl;
    stat_br_ofs<<"\"max_not_solved_branching_rules_bf\":"<<max_not_solved_branching_rules_bf<<","<<endl;
    stat_br_ofs<<"\"not_solved_branching_rules\":"<<not_solved_branching_rules.size()<<","<<endl;
    stat_br_ofs<<"\"time_s\":"<<_time_s<<endl;
    stat_br_ofs<<"}"<<endl;
    stat_br_ofs.close();

    std::ofstream baf_br_ofs (OUTPUT_FILES_PREFIX + ".not_solved_branching_rules.out", std::ofstream::out);
    for(const auto & br : not_solved_branching_rules){
        baf_br_ofs<<export_branching_rule_to_jsonl(br)<<endl;
    }
    baf_br_ofs.close();
    std::ofstream br_ofs (OUTPUT_FILES_PREFIX + ".branching_rules.out", std::ofstream::out);
    for(const auto & br : branching_rules){
        br_ofs<<export_branching_rule_to_jsonl(br)<<endl;
    }
    br_ofs.close();


    proof_close_proof_file();

    return 0;
}