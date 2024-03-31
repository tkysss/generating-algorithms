#pragma once
#include "../branching_rules_computation/branching_rules_primitives.hpp"

std::string branches_to_json(const vector<vector<int>> & branches);

std::string export_branching_rule_to_jsonl(const BranchingRule & br) {
    std::string g_str="";
    g_str+="{\"n\":" + std::to_string(br.g.n) + ",\"edges\":[";
    auto g_edges = br.g.get_edges_list();
    for(int i = 0; i < g_edges.size(); ++i) {
        auto e = g_edges[i];
        if(i!=0)g_str+=",";
        g_str+="["+std::to_string(e.first)+","+std::to_string(e.second)+"]";
    }
    g_str+="]}";

    std::string red_vertices_str="";
    red_vertices_str+="[";
    for(int i = 0; i < br.red_vertices.size(); ++i) {
        if(i!=0)red_vertices_str+=",";
        red_vertices_str+=std::to_string(br.red_vertices[i]);
    }
    red_vertices_str+="]";


    std::string br_str = std::string("{")
        + "\"graph\":" + g_str
        + ",\"red_vertices\":" + red_vertices_str
        + ",\"minimal_branches\":" + branches_to_json(br.minimal_branches)
        + ",\"dominance_free_branches\":" + branches_to_json(br.dominance_free_branches)
        + ",\"branches\":" + branches_to_json(br.adjusted_branches)
        + ",\"bf\":" + std::to_string(br.bf)
        + ",\"type\":" + "\"" + br.type + "\""
    + "}";
    return br_str;
}

std::string branches_to_json(const vector<vector<int>> & branches) {
    std::string branches_str="";
    branches_str+="[";
    for(int i = 0; i < branches.size(); ++i) {
        const auto & branch = branches[i];
        if(i!=0)branches_str+=",";
        std::string branch_str = "";
        branch_str="[";
        for(int j = 0; j < branch.size(); ++j) {
            if(j!=0)branch_str+=",";
            branch_str+=std::to_string(branch[j]);
        }
        branch_str+="]";
        branches_str+=branch_str;
    }
    branches_str+="]";
    return branches_str;
}
