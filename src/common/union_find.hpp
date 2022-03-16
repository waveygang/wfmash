#pragma once

#include <cstdio>
#include <cassert>
#include <vector>
#include <iostream>

namespace union_find {

    // initialize our backing vector
    std::vector<uint64_t> init(uint64_t N) {
        return std::vector<uint64_t>(N, -1);
    }

    // find root of set containing x
    uint64_t root(std::vector<uint64_t>& uf, uint64_t x) {
        uint64_t r = x;
        // find root
        while(uf[r] >= 0)
            r = uf[r];
        // compress path to root
        while(uf[x] >= 0) {
            int tmp = uf[x];
            uf[x] = r;
            x = tmp;
        }
        return r;
    }

    // union of sets containing x and y
    void join(std::vector<uint64_t>& uf, uint64_t x, uint64_t y) {
        x = root(uf, x);
        y = root(uf, y);
        if(x != y) {
            if(uf[x] < uf[y]) {
                uf[x] += uf[y];
                uf[y] = x;
            }
            else {
                uf[y] += uf[x];
                uf[x] = y;
            }
        }
    }

}
