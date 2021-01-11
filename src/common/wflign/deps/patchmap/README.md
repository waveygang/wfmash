# patchmap
A fast and memory efficient hashmap using sorting to resolve collisions

It achieves a very good trade-off between memory efficiency and speed for load factors over ~0.5.
It can be an almost drop-in replacement for std::unordered_map, with the caveat however that
iterators are possibly invalidated by insertions.

Usage:
```C++
#include <iostream>
#include "patchmap.hpp"
int main() {
  whash::patchmap<int,int> hash_table;
  hash_table[7] = 77;
  for (const auto& elem : hash_table) {
    std::cout << elem.first << " " << elem.second << std::endl;
  }
}
```
```bash
> g++ -O3 -std=c++17 -DNDEBUG main.cpp
```

This hashmap is inspired by the 1973 publication "Ordered hash tables".
The idea presented there is to resolve collisions via an ordering defined by the keys directly.
The patchmap however resolves the collisions with an ordering defined by the
hash-value of the keys.
So long there are no collisions the keys in a hash table are already in hash-order.
When there are collisions in the patchmap this order is upheld by inserting
the key-value pairs at the appropriate position, displacing keys that compare greater
to the right and keys that compare less to the left, essentially performing a single
step of insertion sort.
This improves the commonly encountered worst-case complexity of O(n) for lookups to
O(log(n)) as a binary search can always be employed on a ordered list with random access.
In addition, when the hash values are evenly distributed, which is the case for the hash
functions supplied with this library, this allowes for an asymptotic complexity of
O(log(-log(1-a))) because an interpolation search can be used to retrieve the entries,
performing exceptionally well even when the table is almost full.

The resulting ordering is very similar to the one resulting from linear bidirectional
probing with robin hood hashing (see for example sherwood_map) but allowing
interpolation search leads to improved upper bounds while retaining the same average
complexity and overall performance, exceeding other implementations at high load
factors around 96% at the cost of O((1-a)¯²) reordering operations for insertions
and deletions however.

If you are interested using this container contact me and I will make it work for you,
if it does not just work out of the box. No, just contact me regardless, I'm curious
how it works for you!

If you are interested in understanding the patchmap or want to implement it in your
favourite programming language you should have a look at patchmap_v0.hpp.
This is a oversimplified prototype, to ease the understanding without templates
and with binary search instead of interpolation search.

# TODO

 - use boosts advanced allocation to make use of expansion and reallocation
 - re-unify sparse_patchmap.hpp and patchmap.hpp
