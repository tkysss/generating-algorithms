#pragma once

template<int bits>
struct _Bitset {
    static const int N = (bits+63)/64;
    static const int capacity = N*64;
    uint64_t data[N];

    _Bitset() {
        for(int i=0;i<N;++i)data[i]=0;
    }

    void set_bit(int i) {
        assert(i<capacity);
        data[i>>6] |= ((uint64_t)1)<<(i&63);
    }

    void unset_bit(int i) {
        assert(i<capacity);
        data[i>>6] &= ~(((uint64_t)1)<<(i&63));
    }

    bool has_bit(int i) const {
        assert(i<capacity);
        return data[i>>6] & ((uint64_t)1)<<(i&63);
    }

    bool operator < (const _Bitset & o) const {
        for(int i = N - 1; i >= 0; --i) {
            if(data[i]!=o.data[i]) {
                return data[i] < o.data[i];
            }
        }
        return false;
    }

    bool operator == (const _Bitset & o) const {
        for(int i=0;i<N;++i){
            if(data[i]!=o.data[i]) return 0;
        }
        return 1;
    }

    bool operator != (const _Bitset & o)  const {
        return !(*this == o);
    }
};
