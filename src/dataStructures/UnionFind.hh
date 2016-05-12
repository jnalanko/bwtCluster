#ifndef UNIONFIND_HH
#define UNIONFIND_HH

#include <vector>
#include <map>
#include <unordered_map>
#include <stdint.h>

class UnionFind {
    private:
        class UF_element {
            public:
                int32_t size; // group sizes are limited to 80 in Metacluster
                int64_t parent;
                UF_element() : size(1) {}
        };

        std::vector<UF_element> elements;

    public:
        void init(int64_t nElements); // Initializes the structure with elements having ids 0,...,nElements-1
        int64_t find(int64_t id); // Returns the identifier of the set that contains the parameter id
        void doUnion(int64_t id_1, int64_t id_2); // Merges the two sets given by the identifiers
        int64_t getSize(int64_t id); // Returns the size of the set

        UnionFind();
};

#endif
