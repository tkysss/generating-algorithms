#pragma once
#include "headers.hpp"
#include "primitive_operations.hpp"
#include "comparable_bitset.hpp"

extern "C" {
#define MAXN MAX_GRAPH_N
#include "../nauty27rc2/nauty.h"
}

const int _NautyCanonicalRepEdgesBits = MAXN*(MAXN-1)/2;

using _NautyCanonicalRepEdges = _Bitset<_NautyCanonicalRepEdgesBits>;

struct _NautyCanonicalRep {
    int n;
    _NautyCanonicalRepEdges edges;

    explicit _NautyCanonicalRep(): n(0), edges(_NautyCanonicalRepEdges()) {}
    _NautyCanonicalRep(int n, _NautyCanonicalRepEdges edges): n(n), edges(edges) {}

    bool operator < (const _NautyCanonicalRep & o) const {
        if(n != o.n) return n < o.n;
        return edges < o.edges;
    }

    bool operator == (const _NautyCanonicalRep & o) const {
        return n == o.n && edges == o.edges;
    }

    bool operator != (const _NautyCanonicalRep & o) const {
        return !(*this == o);
    }
};

using _NautyRelabeling = std::vector<int>;

struct _NautyCanonicalRepWithRelab {
    _NautyCanonicalRep canonical_rep;
    _NautyRelabeling relabeling;
};

_NautyCanonicalRep make_empty_nauty_canonical_rep() {
    return _NautyCanonicalRep();
}

_NautyCanonicalRepWithRelab make_empty_nauty_canonical_rep_with_relabeling() {
    return {make_empty_nauty_canonical_rep(), {}};
}


// struct _NautyCanonicalRep {
//     int n_size;
//     static const int N = (MAXN*(MAXN-1)/2+63)/64;
//     uint64_t data[N];

//     _NautyCanonicalRep(int n_size) : n_size(n_size) {
//         for(int i = 0; i < N; ++i)data[i]=0;
//     }

//     bool operator < (const _NautyCanonicalRep & o) const {
//         if(n_size != o.n_size) return n_size < o.n_size;
//         for(int i = N - 1; i >= 0; --i) {
//             if(data[i]!=o.data[i]) {
//                 return data[i] < o.data[i];
//             }
//         }
//         return false;
//     }

//     bool operator == (const _NautyCanonicalRep & o) const {
//         if(n_size != o.n_size) return 0;
//         for(int i = 0; i < N; ++i){
//             if(data[i]!=o.data[i]) return 0;
//         }
//         return 1;
//     }

//     bool operator != (const _NautyCanonicalRep & o)  const {
//         return !(*this == o);
//     }

//     void set_bit(int i) {
//         data[i>>6] |= ((uint64_t)1)<<(i&63);
//     }

//     void unset_bit(int i) {
//         data[i>>6] &= ~(((uint64_t)1)<<(i&63));
//     }


//     bool has_bit(int i) const {
//         return data[i>>6] & ((uint64_t)1)<<(i&63);
//     }

//     bool is_valid() const {
//         return n_size != 0;
//     }
// };

const int __nauty_workspace_size = 4*50*MAXM;
struct _Nauty {
    graph __nauty_g[MAXN*MAXM];
    graph __nauty_canong[MAXN*MAXM];
    setword __nauty_workspace[__nauty_workspace_size];
    int __nauty_lab[MAXN],__nauty_ptn[MAXN],__nauty_orbits[MAXN], __nauty_relab[MAXN];
    DEFAULTOPTIONS_GRAPH(__nauty_options);
    statsblk __nauty_stats;
    _Nauty() {
        __nauty_options.getcanon = true;
    }
    _NautyCanonicalRepWithRelab get_canonical_rep(int n, const vector<pair<int,int>> & edges) {
        assert(n<=MAXN);
        int m = SETWORDSNEEDED(n);
        for(int i = 0; i < MAXN*MAXM; ++i)__nauty_canong[i]=0;
        for(int i = 0; i < MAXN; ++i)__nauty_lab[i]=0;
        EMPTYGRAPH(__nauty_g,m,n);
        for(pair<int,int> e : edges)ADDONEEDGE(__nauty_g,e.first,e.second,m);
        nauty(__nauty_g,__nauty_lab,__nauty_ptn,NULL,__nauty_orbits,&__nauty_options,&__nauty_stats,
            __nauty_workspace, __nauty_workspace_size, m,n,__nauty_canong);
        _NautyRelabeling relabeling(n);
        for(int i = 0; i < n; ++i)relabeling[__nauty_lab[i]] = i;
        _NautyCanonicalRepEdges canonical_rep_edges;
        for(pair<int,int> e : edges) {
            int i = relabeling[e.first];
            int j = relabeling[e.second];
            if(i>j)std::swap(i,j);
            canonical_rep_edges.set_bit((j*(j-1)/2+i));
        }
        _NautyCanonicalRep canonical_rep(n, canonical_rep_edges);

        return {canonical_rep, relabeling};
    }
};