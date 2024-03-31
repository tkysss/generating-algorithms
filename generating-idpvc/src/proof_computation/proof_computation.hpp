#pragma once
#include <iostream>
#include <bitset>
#include "../graph_computation/graph_primitives.hpp"
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#ifdef NOPROOF
#define PROOF_CONTROL return;
#else
#define PROOF_CONTROL
#endif



boost::iostreams::filtering_streambuf<boost::iostreams::output> * proof_outbuf = 0;
std::ofstream * proof_out_file;
std::ostream proof_out(0);

void proof_open_proof_file() {
    PROOF_CONTROL;
    proof_out_file = new std::ofstream(OUTPUT_FILES_PREFIX+".proof.out.gz", std::ios_base::out | std::ios_base::binary);

    proof_outbuf = new boost::iostreams::filtering_streambuf<boost::iostreams::output>();
    proof_outbuf->push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(9)));
    proof_outbuf->push(*proof_out_file);
    proof_out.rdbuf(proof_outbuf);
}

void proof_close_proof_file() {
    PROOF_CONTROL;
    delete proof_outbuf;
    proof_out_file->close();
}

void proof_parameters(int DPVC_PATH_LEN, double DPVC_BF, int GENERATIONS) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"parameters\"";
    proof_out<<",";
    proof_out<<"\"dpvc_path_len\":"<<DPVC_PATH_LEN;
    proof_out<<",";
    proof_out<<"\"dpvc_bf\":"<<DPVC_BF;
    proof_out<<",";
    proof_out<<"\"generations\":"<<GENERATIONS;
    proof_out<<"}";
    proof_out<<endl;
}

void _proof_output_vector_json(const vector<int> & v, bool convert_int_to_vector=false) {
    PROOF_CONTROL;
    proof_out<<"[";
    bool not_first = false;
    for(int e : v) {
        if(not_first)proof_out<<",";
        not_first=true;
        if(convert_int_to_vector) {
            _proof_output_vector_json(bit_mask_to_vector(e), false);
        }
        else {
            proof_out<<e;
        }
    }
    proof_out<<"]";
}

void _proof_output_vector_of_vector_json(const vector<vector<int>> & vv, bool convert_int_to_vector=false) {
    PROOF_CONTROL;
    proof_out<<"[";
    bool not_first = false;
    for(const vector<int> v : vv) {
        if(not_first)proof_out<<",";
        not_first=true;
        _proof_output_vector_json(v, convert_int_to_vector);
    }
    proof_out<<"]";
}

template<int _Bits>
void _proof_output_bitset_json(const _Bitset<_Bits> & _bitset) {
    PROOF_CONTROL;
    proof_out<<"\"";
    std::string bitset_out;
    bool got_non_zero=false;
    for(int i = _bitset.N - 1; i >= 0; --i) {
        if(_bitset.data[i] != 0) {
            got_non_zero=true;
        }
        if(got_non_zero) {
            bitset_out += std::bitset<64>(_bitset.data[i]).to_string();
        }
    }
    auto first_one = bitset_out.find_first_of('1');
    if(first_one == std::string::npos) {
        // Should not happen, no empty graphs allowed.
        assert(0);
    }
    else {
        bitset_out = bitset_out.substr(first_one);
    }

    proof_out<<bitset_out;
    proof_out<<"\"";
}

void _proof_output_canonical_rep_with_relabeling_json(const _NautyCanonicalRepWithRelab & canonical_rep_with_relabeling) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"c_n\":"<<canonical_rep_with_relabeling.canonical_rep.n;
    proof_out<<",";
    proof_out<<"\"c_edges\":";_proof_output_bitset_json(canonical_rep_with_relabeling.canonical_rep.edges);
    proof_out<<",";
    proof_out<<"\"c_relab\":";_proof_output_vector_json(canonical_rep_with_relabeling.relabeling);
    proof_out<<"}";
}

void _proof_output_graph_json(const _Graph & g) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"n\":"<<g.n;
    proof_out<<",";
    proof_out<<"\"edges\":";_proof_output_bitset_json(g.edges);
    proof_out<<",";
    proof_out<<"\"canon\":";_proof_output_canonical_rep_with_relabeling_json(g.canonical_rep_with_relabeling);
    proof_out<<",";
    proof_out<<"\"red_vs\":";_proof_output_vector_json(bit_mask_to_vector(g.red_vertices_mask));

    proof_out<<"}";
}

void proof_current_bumpy_graphs(const vector<_Graph> & current_bumpy_graphs) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"current_bumpy_graphs_start\"";
    proof_out<<"}";
    proof_out<<endl;

    for(const _Graph & g : current_bumpy_graphs) {
        proof_out<<"{";
        proof_out<<"\"t\":\"current_bumpy_graph\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<"}";
        proof_out<<endl;
    }

    proof_out<<"{";
    proof_out<<"\"t\":\"current_bumpy_graphs_end\"";
    proof_out<<"}";
    proof_out<<endl;
}

void proof_current_bumpy_graphs_with_expansions(const vector<_Graph> & current_bumpy_graphs, const vector<vector<_Graph>> & current_bumpy_graph_expansions) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"current_bumpy_graphs_with_expansions_start\"";
    proof_out<<"}";
    proof_out<<endl;

    for(int i = 0; i < current_bumpy_graphs.size(); ++i) {
        const _Graph & g = current_bumpy_graphs[i];
        proof_out<<"{";
        proof_out<<"\"t\":\"current_bumpy_graph\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<"}";
        proof_out<<endl;

        proof_out<<"{";
        proof_out<<"\"t\":\"current_bumpy_graph_expansions_start\"";
        proof_out<<"}";
        proof_out<<endl;

        for(const _Graph & eg : current_bumpy_graph_expansions[i]) {
            proof_out<<"{";
            proof_out<<"\"t\":\"current_bumpy_graph_expansion\"";
            proof_out<<",";
            proof_out<<"\"eg\":";_proof_output_graph_json(eg);
            proof_out<<"}";
            proof_out<<endl;
        }

        proof_out<<"{";
        proof_out<<"\"t\":\"current_bumpy_graph_expansions_end\"";
        proof_out<<"}";
        proof_out<<endl;
    }

    proof_out<<"{";
    proof_out<<"\"t\":\"current_bumpy_graphs_with_expansions_end\"";
    proof_out<<"}";
    proof_out<<endl;
}

void proof_generation_start(int generation, const vector<_Graph> & current_bumpy_graphs) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"generation_start\"";
    proof_out<<",";
    proof_out<<"\"generation\":"<<generation;
    proof_out<<"}";
    proof_out<<endl;
    proof_current_bumpy_graphs(current_bumpy_graphs);
}

void proof_generation_end(int generation, const vector<_Graph> & current_bumpy_graphs, const vector<vector<_Graph>> & current_bumpy_graph_expansions) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"generation_end\"";
    proof_out<<",";
    proof_out<<"\"generation\":"<<generation;
    proof_out<<"}";
    proof_out<<endl;
    proof_current_bumpy_graphs_with_expansions(current_bumpy_graphs, current_bumpy_graph_expansions);
}

void proof_generation_next(int generation, const vector<_Graph> & next_bumpy_graphs) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"generation_next\"";
    proof_out<<",";
    proof_out<<"\"generation\":"<<generation;
    proof_out<<"}";
    proof_out<<endl;
    proof_current_bumpy_graphs(next_bumpy_graphs);
}

void proof_expansion_start() {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"expansion_start\"";
    proof_out<<"}";
    proof_out<<endl;
}

void proof_expansion_end() {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"expansion_end\"";
    proof_out<<"}";
    proof_out<<endl;
}

void proof_round_start(int round) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"round_start\"";
    proof_out<<",";
    proof_out<<"\"round\":"<<round;
    proof_out<<"}";
    proof_out<<endl;
}

void proof_round_end(int round) {
    PROOF_CONTROL;
    proof_out<<"{";
    proof_out<<"\"t\":\"round_end\"";
    proof_out<<",";
    proof_out<<"\"round\":"<<round;
    proof_out<<"}";
    proof_out<<endl;
}

void proof_expand_graph_and_filter_forbidden_subgraphs_filtered(const _Graph & g, const _Graph & eg, int _proof_match_vertices_mask, const _NautyCanonicalRepWithRelab & _proof_match_canonical_rep_with_relabeling) {
    PROOF_CONTROL;
    #pragma omp critical
    {
    if(_proof_match_vertices_mask == -1) {
        proof_out<<"{";
        proof_out<<"\"t\":\"expand_gaffsf_cache\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<",";
        proof_out<<"\"eg\":";_proof_output_graph_json(eg);
        proof_out<<"}";
        proof_out<<endl;
    }
    else {
        proof_out<<"{";
        proof_out<<"\"t\":\"expand_gaffsf_match\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<",";
        proof_out<<"\"eg\":";_proof_output_graph_json(eg);
        proof_out<<",";
        proof_out<<"\"m_vs\":";_proof_output_vector_json(bit_mask_to_vector(_proof_match_vertices_mask));
        proof_out<<",";
        proof_out<<"\"m_canon\":";_proof_output_canonical_rep_with_relabeling_json(_proof_match_canonical_rep_with_relabeling);
        proof_out<<"}";
        proof_out<<endl;
    }
    }
}

void proof_handled_by_red_component_reduction(const _Graph & g, const vector<int> & red_component1, const vector<int> & red_component2, int red_component_star_vertex) {
    PROOF_CONTROL;
    #pragma omp critical
    {
        proof_out<<"{";
        proof_out<<"\"t\":\"handled_by_red_component_reduction\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<",";
        proof_out<<"\"rc1\":";_proof_output_vector_json(red_component1);
        proof_out<<",";
        proof_out<<"\"rc2\":";_proof_output_vector_json(red_component2);
        proof_out<<",";
        proof_out<<"\"v\":"<<red_component_star_vertex;
        proof_out<<"}";
        proof_out<<endl;
    }
}

void proof_handled_by_red_star_reduction(const _Graph & g, int kcenter_vertices_mask, int red_rays_vertices_mask) {
    PROOF_CONTROL;
    #pragma omp critical
    {
        proof_out<<"{";
        proof_out<<"\"t\":\"handled_by_red_star_reduction\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<",";
        proof_out<<"\"kcenter_vertices\":";_proof_output_vector_json(bit_mask_to_vector(kcenter_vertices_mask));
        proof_out<<",";
        proof_out<<"\"red_rays_vertices\":";_proof_output_vector_json(bit_mask_to_vector(red_rays_vertices_mask));
        proof_out<<"}";
        proof_out<<endl;
    }
}

void proof_handled_by_vertex_cover_struction(const _Graph & g, int v, int p, int q) {
    PROOF_CONTROL;
    #pragma omp critical
    {
        proof_out<<"{";
        proof_out<<"\"t\":\"handled_by_vertex_cover_struction\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<",";
        proof_out<<"\"v\":"<<v;
        proof_out<<",";
        proof_out<<"\"p\":"<<p;
        proof_out<<",";
        proof_out<<"\"q\":"<<q;
        proof_out<<"}";
        proof_out<<endl;
    }
}

void proof_handled_by_vertex_cover_dominance(const _Graph & g, int u, int v) {
    PROOF_CONTROL;
    #pragma omp critical
    {
        proof_out<<"{";
        proof_out<<"\"t\":\"handled_by_vertex_cover_dominance\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<",";
        proof_out<<"\"u\":"<<u;
        proof_out<<",";
        proof_out<<"\"v\":"<<v;
        proof_out<<"}";
        proof_out<<endl;
    }
}

void proof_branching_rule(const _Graph & g, const vector<vector<int>> & minimal_solutions, const vector<vector<int>> & dominance_free_solutions, const pair<vector<vector<int>>, vector<vector<int>>> & solution_dominance_graph_adj_r, const vector<vector<int>> & adjusted_solutions) {
    PROOF_CONTROL;
    #pragma omp critical
    {
        proof_out<<"{";
        proof_out<<"\"t\":\"branching_rule\"";
        proof_out<<",";
        proof_out<<"\"g\":";_proof_output_graph_json(g);
        proof_out<<",";
        proof_out<<"\"minimal_solutions\":";_proof_output_vector_of_vector_json(minimal_solutions);
        proof_out<<",";
        proof_out<<"\"dominance_free_solutions\":";_proof_output_vector_of_vector_json(dominance_free_solutions);
        proof_out<<",";
        proof_out<<"\"solution_dominance_graph_adj\":";_proof_output_vector_of_vector_json(solution_dominance_graph_adj_r.first);
        proof_out<<",";
        proof_out<<"\"solution_dominance_graph_adj_red_vertices\":";_proof_output_vector_of_vector_json(solution_dominance_graph_adj_r.second, true);
        proof_out<<",";
        proof_out<<"\"adjusted_solutions\":";_proof_output_vector_of_vector_json(adjusted_solutions);
        proof_out<<"}";
        proof_out<<endl;
    }
}
