#pragma once
#include "headers.hpp"

vector<vector<int>> ordered_combinations(const vector<int> & src, int k) {
    assert(0<=k);
    assert(k<=src.size());
    vector<vector<int>> out;
    vector<int> idx(k);
    for(int i = 0; i < k; ++i) {
        idx[i]=i;
    }
    while(1){
        vector<int> x;
        for(int i = 0; i < k; ++i) {
            x.push_back(src[idx[i]]);
        }
        out.push_back(move(x));
        int i;
        for(i=k-1;i>=0;--i) {
            if(idx[i] < src.size() - k + i) {
                ++idx[i];
                for(int j=i+1; j < k; ++j)idx[j]=idx[j-1]+1;
                break;
            }
        }
        if(i<0)break;
    }
    return out;
}

vector<int> ordered_combinations_of_bits(int n, int k) {
    assert(n<=31);
    assert(0<=k);
    assert(k<=n);
    if(k==0){return {0};}
    const int nmask = (1<<n)-1;
    vector<int> out;
    int c = (1<<k)-1;
    while(1){
        out.push_back(c);
        int idx_of01 = __builtin_ctz(c);
        int cnt_of1 = 0;
        for(; idx_of01 < n; ++idx_of01) {
            int m = (1<<idx_of01) + (1<<(idx_of01+1));
            if((m&c) != m)break;
            cnt_of1++;
        }
        int uppermask = nmask&~((1<<(idx_of01+2))-1);
        c = (c&uppermask) | (1<<(idx_of01+1)) | ((1<<cnt_of1)-1);
        if(c&(1<<n))break;
    }
    return out;
}

vector<vector<int>> ordered_powerset_lo_hi(const vector<int> & src, int lo, int hi) {
    assert(0<=lo);
    assert(lo<=hi);
    assert(hi<=src.size());
    vector<vector<int>> out;
    for(int i=lo;i<=hi;++i)
        for(vector<int> & c : ordered_combinations(src, i))
            out.push_back(move(c));
    return out;
}

vector<vector<int>> ordered_powerset(const vector<int> & src) {
    return ordered_powerset_lo_hi(src, 0, src.size());
}

vector<vector<int>> ordered_powerset_nonempty(const vector<int> & src) {
    return ordered_powerset_lo_hi(src, 1, src.size());
}

vector<vector<int>> ordered_powerset_nonfull(const vector<int> & src) {
    return ordered_powerset_lo_hi(src, 0, src.size()-1);
}

vector<int> range(int n) {
    assert(n>=1);
    vector<int> x(n);for(int i = 0; i < n; ++i)x[i]=i;return x;
}

vector<int> set_union(const vector<int> & a, const vector<int> & b) {
    vector<int> r(b);
    for(int aa : a){
        bool ok=1;
        for(int bb : b) {
            if(aa==bb){ok=0;break;}
        }
        if(ok)r.push_back(aa);
    }
    sort(r.begin(), r.end());
    return r;
}

vector<int> set_minus(const vector<int> & a, const vector<int> & b) {
    vector<int> r;
    for(int aa : a){
        bool ok=1;
        for(int bb : b) {
            if(aa==bb){ok=0;break;}
        }
        if(ok)r.push_back(aa);
    }
    sort(r.begin(), r.end());
    return r;
}

vector<int> set_intersection(const vector<int> & a, const vector<int> & b) {
    vector<int> r;
    for(int aa : a){
        bool ok=0;
        for(int bb : b) {
            if(aa==bb){ok=1;break;}
        }
        if(ok)r.push_back(aa);
    }
    sort(r.begin(), r.end());
    return r;
}

bool is_in_set(int x, const vector<int> & s) {
    for(int ss : s)if(x==ss)return true;
    return false;
}

bool is_subset(const vector<int> & a, const vector<int> & b) {
    return set_minus(a,b).size() == 0;
}

vector<int> sorted(const vector<int> & a) {
    vector<int> r(a);
    sort(r.begin(), r.end());
    return r;
}

int vector_to_bit_mask(const vector<int> & vi) {
    int bit_mask = 0;
    for(int v : vi) {
        assert(0<=v);
        assert(v<=30);
        bit_mask |= 1<<v;
    }
    return bit_mask;
}


vector<int> convert_vectors_to_bit_masks(const vector<vector<int>> & vvi) {
    vector<int> bit_masks;
    for(const vector<int> & vi : vvi) {
        bit_masks.push_back(vector_to_bit_mask(vi));
    }
    return bit_masks;
}

int size_bit_masks(int a) {
    return __builtin_popcount(a);
}

int set_intersection_bit_masks(int a, int b) {
    return a&b;
}

int set_union_bit_masks(int a, int b) {
    return a|b;
}

bool is_subset_bit_masks(int a, int b) {
    return (a&b)==a;
}

int set_minus_bit_masks(int a, int b) {
    return a&(~b);
}

vector<int> bit_mask_to_vector(int a) {
    assert(a>=0);
    vector<int> r;
    int i = 0;
    while(a) {
        if(a&1)r.push_back(i);
        ++i;
        a>>=1;
    }
    return r;
}

vector<vector<int>> convert_bit_masks_to_vectors(const vector<int> & vi) {
    vector<vector<int>> r;
    for(int b : vi) r.push_back(bit_mask_to_vector(b));
    return r;
}
