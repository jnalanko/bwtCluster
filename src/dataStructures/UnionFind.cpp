#include "UnionFind.hh"
#include <utility>
#include <iostream>

using namespace std;

UnionFind::UnionFind(){}

void UnionFind::init(int64_t nElements){
    elements.resize(nElements); elements.shrink_to_fit();
    for(int64_t i = 0; i < nElements; i++) elements[i].parent = i;
}

int64_t UnionFind::find(int64_t id){
    if(id < 0) cout << "Negative id: " << id << endl;
    int64_t x = id;
    while(x != elements[x].parent)
        x = elements[x].parent;
    return x;
}

void UnionFind::doUnion(int64_t id_1, int64_t id_2){
    int64_t p1 = find(id_1);
    int64_t p2 = find(id_2);
    if(p1 == p2) return;
    if(elements[p1].size < elements[p2].size){
        elements[p1].parent = p2;
        elements[p2].size += elements[p1].size;
    }
    else if(elements[p2].size < elements[p1].size){
        elements[p2].parent = p1;
        elements[p1].size += elements[p2].size;
    }
    else{ // Equal sizes
        elements[p2].parent = p1;
        elements[p1].size += elements[p2].size;
    }
}

int64_t UnionFind::getSize(int64_t id){
    return elements[find(id)].size;
}
