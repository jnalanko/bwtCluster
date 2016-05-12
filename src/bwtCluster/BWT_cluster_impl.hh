#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <ctime>
#include <stdexcept>
#include "BWT_cluster.hh"
#include "Parser.hh"
#include "tools.hh"
#include <cassert>
#include "UnionFind.hh"
#include <sdsl/bit_vectors.hpp>
#include "BD_BWT_index.hh"
#include <chrono>
#include <set>
#include "malloc.h"
#include "Read_group_builder.hh"
    
using namespace std;

template<class t_bitvector> 
BWT_cluster<t_bitvector>::BWT_cluster(string workspace_path, std::string paired_end_file, int64_t nThreads) 
    : nReads(-1), space_bits(0), start_time(getUnixEpoch()), nThreads(nThreads), workspace_path(workspace_path), paired_end_file(paired_end_file),
    memory_log(workspace_path + "/memory.log") {}
   
template<class t_bitvector> 
BWT_cluster<t_bitvector>::~BWT_cluster() {}

template<class t_bitvector> 
void BWT_cluster<t_bitvector>::memlog(int64_t space, string label){
    memory_log << getUnixEpoch() - start_time << " " << space << " " << label << endl;
}

// Find the bwt position of the first suffix of the given read (not the dollar)
template<class t_bitvector> 
int64_t BWT_cluster<t_bitvector>::locate_read_start(const string& read){

    Interval_pair I(0,dataLength-1,0,dataLength-1);
    for(int64_t i = read.size() - 1; i >= 0; i--){ 
        I = left_extend(index, I, read[i]);
    }
    I = left_extend(index, I, '$');
    if(I.forward.size() > 1){
        cerr << "Error: duplicate read in data: " << read << endl;
        exit(1);
    }
    assert(I.forward.size() == 1);
    assert(I.forward.left == I.forward.right);
    return I.forward.left;
}

// NOTE: nReads must be set before calling this
// NOTE: index must be built before calling this
template<class t_bitvector>
std::vector<BWT_inversion_iterator<t_bitvector> > BWT_cluster<t_bitvector>::get_iterators_sorted(){
    assert(nReads >= nThreads);
    std::vector<BWT_inversion_iterator<t_bitvector> > iterators;
    vector<int64_t> read_ids;
    for(int64_t i = 0; i < nThreads; i++){
        read_ids.push_back((nReads / nThreads) * i); 
    }

    for(int64_t i = 0; i < nThreads; i++){
        // The block of dollars in the BWT contains the read start suffixes
        // in sorted order. The position of read i is i + 1.
        int64_t first_pos = (i == nThreads - 1) ? 0 : read_ids[(i + 1) % nThreads] + 1;
        int64_t last_pos = read_ids[i] + 1;
        int64_t first_id = read_ids[(i + 1) % nThreads] - 1;
        if(first_id == -1) first_id = nReads - 1;
        BWT_inversion_iterator<t_bitvector> it(&index.forward_bwt, first_pos, last_pos, first_id);
        iterators.push_back(it);
    }
    return iterators;
}

template<class t_bitvector>
void BWT_cluster<t_bitvector>::loadData(std::string fasta_filename){
    // Clear old index
    BD_BWT_index<t_bitvector> new_index;
    this->index = new_index;
    
    cerr << getTimeString() << " Parsing input data" << endl;
    // Parse the data
    ifstream input(fasta_filename);
    Parser parser;
    pair<string,int64_t> parseResult = parser.parseConcatenate(input,'$');
    string& concatenation = parseResult.first;
    dataLength = concatenation.size() + 1; // +1 = final dollar appended to get the BWT
    nReads = parseResult.second;
    concatenation = ""; concatenation.shrink_to_fit();
    
    // Build the index
    cerr << getTimeString() << " Initializing the index" << endl;
    buildBWT_ropebwt(fasta_filename);

    cerr << getTimeString() << " Number of characters in data: " << dataLength << endl;
    cerr << getTimeString() << " Number of reads in data: " << nReads << endl;
    
    // Compute and print the size of the index in bytes
    int64_t index_bytes = size_in_bytes(index.forward_bwt) +
                          size_in_bytes(index.reverse_bwt);
    cerr << getTimeString() <<  " Index takes: " << index_bytes << " bytes (" <<
        index_bytes*8 << " bits)" << endl;
        
}

template<class t_bitvector>
vector<int64_t> BWT_cluster<t_bitvector>::filterReads(int k, int tau, bool rc, std::string good_reads, std::string bad_reads){
    
    Filterer<t_bitvector> filterer(&index, get_iterators_sorted(), workspace_path, paired_end_file);
    if(rc)
        return filterer.filterReadsRC(k,tau,good_reads,bad_reads);
    else
        return filterer.filterReads(k,tau,good_reads,bad_reads);

}

template<class t_bitvector>
void BWT_cluster<t_bitvector>::groupReads(int64_t k, bool rc, int64_t maxGroupSize, std::ostream &output){

    Read_group_builder<t_bitvector> rgb(&index, get_iterators_sorted(), paired_end_file);
    
    if(rc)
        rgb.groupReadsRC(k,maxGroupSize,workspace_path,output);
    else
        rgb.groupReads(k,maxGroupSize,output);
    
}

template<class t_bitvector>
void BWT_cluster<t_bitvector>::groupReadsSubmaximal(int64_t k, bool rc, int64_t maxGroupSize, std::ostream &output){
    if(rc){
        cerr << "Error: Reverse-complement submaximal intervals not implemented" << endl;
        exit(1);
    }    
    
    Read_group_builder<t_bitvector> rgb(&index,get_iterators_sorted(), paired_end_file);
    rgb.groupReadsSubmaximal(k,maxGroupSize,output);
}



template<class t_bitvector>
void BWT_cluster<t_bitvector>::filterGroups(std::istream &in, int64_t minGroupSize, std::ostream &newGroups){
    cerr << getTimeString() << " Filtering groups" << endl;
    vector<int64_t> currentGroup;
    string line;
    int64_t reads_kept = 0;
    int64_t number_of_preclusters = 0;
    int64_t number_of_large_preclusters = 0;
    while(getline(in, line)){
        number_of_preclusters++;
        stringstream ss(line);
        int64_t readId;
        while(ss >> readId)
            currentGroup.push_back(readId);
        if(currentGroup.size() >= minGroupSize){
            reads_kept += currentGroup.size();
            number_of_large_preclusters++;
            for(int64_t x : currentGroup){
                newGroups << x << " ";
            }
            newGroups << "\n";
        }
        currentGroup.clear();
    }
    cerr << getTimeString() << " " << number_of_preclusters << " preclusters found" << endl;
    cerr << getTimeString() << " " << number_of_large_preclusters << " preclusters of size >= " << minGroupSize << " found" << endl;
    cerr << getTimeString() << " Kept " << reads_kept << " reads" << endl;

}

template<class t_bitvector>
void BWT_cluster<t_bitvector>::buildBWT_psascan(string& input, string psascan_path){
    if(input.size() == 0){
        throw runtime_error(string("Error at "  + string(__FILE__) + " line " + to_string(__LINE__) +
                    ": Tried to build a BWT of an empty string"));
    }
    
    initialize_index_external_memory_psascan(index,input,workspace_path, psascan_path);
    cerr << getTimeString() << " Index alphabet: ";
    for(char c : index.alphabet) cerr << "'" << c << "' ";
    cerr << endl;
    
}

template<class t_bitvector>
void BWT_cluster<t_bitvector>::buildBWT_ropebwt(string fasta_filename){
    
    initialize_index_external_memory_ropebwt(index,fasta_filename,workspace_path,dataLength/2);
    cerr << getTimeString() << " Index alphabet: ";
    for(char c : index.alphabet) cerr << "'" << c << "' ";
    cerr << endl;
    
}