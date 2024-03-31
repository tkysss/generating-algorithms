#pragma once
#include "../common_primitives/nauty.hpp"
#include <bitset>

//
// Edge storage structure
// 0   1   3   6
// 0,1 0,2 0,3 0,4
// 2   4   7
// 1,2 1,3 1,4
// 5   8
// 2,3 2,4
// 9
// 3,4
//
// We store edges efficiently as a bitset.
const int _GraphEdgesBits = MAX_GRAPH_N*(MAX_GRAPH_N-1)/2;
using _GraphEdges = _Bitset<_GraphEdgesBits>;

struct _Graph {
    int n;
    _GraphEdges edges;
    _NautyCanonicalRepWithRelab canonical_rep_with_relabeling;
    int red_vertices_mask;

    _Graph(int n, const vector<pair<int,int>> & edges, const vector<int> & red_vertices, _Nauty & nt)
        : n(n), canonical_rep_with_relabeling(nt.get_canonical_rep(n, edges)) {
        assert(n<=MAX_GRAPH_N);
        for(auto e : edges)set_edge(e.first, e.second);
        red_vertices_mask = 0;
        for(int rv : red_vertices) red_vertices_mask |= 1<<rv;
    }

    _Graph(int n, const vector<pair<int,int>> & edges, int red_vertices_mask, _Nauty & nt)
        : n(n), canonical_rep_with_relabeling(nt.get_canonical_rep(n, edges)), red_vertices_mask(red_vertices_mask) {
        assert(n<=MAX_GRAPH_N);
        for(auto e : edges)set_edge(e.first, e.second);
    }

    _Graph(int n, _GraphEdges edges, _NautyCanonicalRepWithRelab canonical_rep_with_relabeling, int red_vertices_mask)
        : n(n), edges(edges), canonical_rep_with_relabeling(canonical_rep_with_relabeling), red_vertices_mask(red_vertices_mask) {}

    bool has_edge(int i, int j) const {
        assert(i<j);
        return edges.has_bit((j*(j-1)/2+i));
    }

    void set_edge(int i, int j) {
        assert(i<j);
        edges.set_bit((j*(j-1)/2+i));
    }

    int get_adjacent_vertices_mask(int v) const {
        int av = 0;
        for(int i = 0; i < v; ++i)if(has_edge(i,v))av|=1<<i;
        for(int i=v+1; i < n; ++i)if(has_edge(v,i))av|=1<<i;
        return av;
    }

    vector<int> get_adjacent_vertices(int v) const {
        vector<int> av;
        for(int i = 0; i < v; ++i)if(has_edge(i,v))av.push_back(i);
        for(int i=v+1; i < n; ++i)if(has_edge(v,i))av.push_back(i);
        return av;
    }


    vector<vector<int>> get_adjacent_lists() const {
        vector<vector<int>>  avs(n);
        for(int i=0;i<n;++i){
            for(int j=i+1;j<n;++j){
                if(has_edge(i,j)){
                    avs[i].push_back(j);
                    avs[j].push_back(i);
                }
            }
        }
        return avs;
    }

    bool annotate_red_vertices(const vector<_Graph> & g_expansions) {
        // The new vertex in the expansions for the graph is always "n".
        // Therefore, the vertex v is red only if it is not adjacent to vertex "n" in any of its current g_expansions.

        // To track if something changed.
        int old_red_vertices_mask = red_vertices_mask;

        for(int v = 0; v < n; ++v) {
            // OPTIMIZATION:
            // If the vertex is already red, we do not need to check it again.
            if(red_vertices_mask & (1<<v)) continue;
            bool is_red = true;
            for(const _Graph & eg : g_expansions) {
                if(eg.has_edge(v,n)){is_red = false; break;}
            }
            if(is_red)red_vertices_mask |= 1<<v;
        }

        return old_red_vertices_mask != red_vertices_mask;
    }

    bool operator < (const _Graph & og) const {
        return canonical_rep_with_relabeling.canonical_rep < og.canonical_rep_with_relabeling.canonical_rep;
    }

    vector<pair<int,int>> get_edges_list() const {
        vector<pair<int,int>> el;
        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                if(has_edge(i,j)){el.push_back({i, j});}
            }
        }
        return el;
    }

    _GraphEdges get_induced_edges_mask(int vs_mask) const {
        _GraphEdges ie;
        vector<int> vs_remap(n, -1);
        int vs_cnt=0;
        for(int i = 0; i<n; ++i) {
            if(vs_mask&(1<<i)) {
                vs_remap[vs_cnt] = i;
                vs_cnt++;
            }
        }

        for(int u = 0; u < vs_cnt; ++u) {
            for(int v = u+1; v < vs_cnt; ++v) {
                if(has_edge(vs_remap[u], vs_remap[v]))ie.set_bit((v*(v-1)/2+u));
            }
        }

        return ie;
    }

    vector<pair<int,int>> get_induced_edges_list(int vs_mask) const {
        vector<pair<int,int>> el;
        vector<int> vs_remap(n, -1);
        int vs_cnt=0;
        for(int i = 0; i<n; ++i) {
            if(vs_mask&(1<<i)) {
                vs_remap[vs_cnt] = i;
                vs_cnt++;
            }
        }

        for(int u = 0; u < vs_cnt; ++u) {
            for(int v = u+1; v < vs_cnt; ++v) {
                if(has_edge(vs_remap[u], vs_remap[v]))el.push_back({u, v});
            }
        }

        return el;
    }

    std::string get_json() const {
        std::string g_str="";
        g_str+="{\"n\":" + std::to_string(n) + ",\"edges\":[";
        auto es = get_edges_list();
        for(int i = 0; i < es.size(); ++i) {
            if(i!=0)g_str+=",";
            g_str+="["+std::to_string(es[i].first)+","+std::to_string(es[i].second)+"]";
        }
        g_str+="]}";
        return g_str;
    }
};
