#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <set>
#include <deque>
#include <stack>
#include <map>
#include "bwt.hh"
#include "tools.hh"
#include "malloc.h"
#include "fasta_tools.hh"

using namespace sdsl;
using namespace std;

template<class t_bitvector>
int64_t compute_cumulative_char_rank_in_interval(const wt_huff<t_bitvector>& wt, vector<uint8_t> alphabet, char c, Interval I){
    int64_t ans = 0;
    if(I.size() == 0) return 0;
    
    // Sum of ranks of all characters that are lexicographically smaller than c
    for(char d : alphabet){
        if(d == c) break;
        ans += wt.rank(I.right + 1,d) - wt.rank(I.left,d);
    }
    return ans;
}


template<class t_bitvector>
vector<uint8_t> get_interval_symbols(const wt_huff<t_bitvector>& wt, Interval I){
    if(I.size() == 0){
        vector<uint8_t> empty;
        return empty;
    }
    
    int_vector_size_type nExtensions;
    vector<uint8_t> symbols(wt.sigma);
    vector<uint64_t> ranks_i(wt.sigma);
    vector<uint64_t> ranks_j(wt.sigma);
    wt.interval_symbols(I.left, I.right+1, nExtensions, symbols, ranks_i, ranks_j);
    while(symbols.size() > nExtensions) symbols.pop_back();
    return symbols;
}

template<class t_bitvector>
void get_interval_symbols(const wt_huff<t_bitvector>& wt, Interval I, int_vector_size_type& nExtensions, vector<uint8_t>& symbols,
 vector<uint64_t>& ranks_i, vector<uint64_t>& ranks_j){
    if(I.size() == 0){
        nExtensions = 0;
        return;
    }
    
    wt.interval_symbols(I.left, I.right+1, nExtensions, symbols, ranks_i, ranks_j);
}


// Takes a backward step in the forward bwt
template<class t_bitvector>
int64_t backward_step(const BD_BWT_index<t_bitvector>& index, int64_t pos){
    char c = index.forward_bwt[pos];
    return index.forward_char_block_starts[c] + index.forward_bwt.rank(pos, c); // TODO: untested
}

template<class t_bitvector>
Interval_pair left_extend(const BD_BWT_index<t_bitvector>& index, Interval_pair intervals, char c){
    int64_t cumul_rank_c = compute_cumulative_char_rank_in_interval(index.forward_bwt, index.alphabet, c, intervals.forward);
    return left_extend(index,intervals,c,cumul_rank_c);
}

template<class t_bitvector>
Interval_pair right_extend(const BD_BWT_index<t_bitvector>& index, Interval_pair intervals, char c){
    int64_t cumul_rank_c = compute_cumulative_char_rank_in_interval(index.reverse_bwt, index.alphabet, c, intervals.reverse);
    return right_extend(index,intervals,c,cumul_rank_c);
}

template<class t_bitvector>
Interval_pair left_extend(const BD_BWT_index<t_bitvector>& index, Interval_pair intervals, char c, int64_t cumul_rank_c){
    if(intervals.forward.size() == 0)
        return Interval_pair(-1,-2,-1,-2);
    
    Interval forward = intervals.forward;
    Interval reverse = intervals.reverse;
    
    // Compute the new forward interval
    int64_t num_c_in_interval = index.forward_bwt.rank(forward.right + 1,c) - index.forward_bwt.rank(forward.left,c);
    int64_t start_f_new = index.cumulative_char_count[c] + index.forward_bwt.rank(forward.left, c); // Start in forward
    int64_t end_f_new = start_f_new + num_c_in_interval - 1; // End in forward

    if(start_f_new > end_f_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
    
    // Compute the new reverse interval
    int64_t start_r_new = reverse.left + cumul_rank_c;
    int64_t end_r_new = start_r_new + (end_f_new - start_f_new); // The forward and reverse intervals must have same length
    
    return Interval_pair(start_f_new,end_f_new,start_r_new,end_r_new);
}

template<class t_bitvector>
Interval_pair right_extend(const BD_BWT_index<t_bitvector>& index, Interval_pair intervals, char c, int64_t cumul_rank_c){
    if(intervals.forward.size() == 0)
        return Interval_pair(-1,-2,-1,-2);
    
    Interval forward = intervals.forward;
    Interval reverse = intervals.reverse;
    
    // Compute the new reverse interval
    int64_t num_c_in_interval = index.reverse_bwt.rank(reverse.right + 1,c) - index.reverse_bwt.rank(reverse.left,c);
    int64_t start_r_new = index.cumulative_char_count[c] + index.reverse_bwt.rank(reverse.left, c); // Start in reverse
    int64_t end_r_new = start_r_new + num_c_in_interval - 1; // End in reverse

    if(start_r_new > end_r_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
    
    // Compute the new forward interval
    int64_t cumul_rank = cumul_rank_c;
    int64_t start_f_new = forward.left + cumul_rank;
    int64_t end_f_new = start_f_new + (end_r_new - start_r_new); // The forward and reverse intervals must have same length
    
    return Interval_pair(start_f_new,end_f_new,start_r_new,end_r_new);
    
}


// Compute the cumulative sum of character counts in lexicographical order
// Assumes alphabet is sorted
// Counts = vector with 256 elements
template<class t_bitvector>
void count_smaller_chars(const wt_huff<t_bitvector>& bwt, vector<uint8_t>& alphabet, vector<int64_t>& counts, Interval I){
    assert(alphabet.size() != 0);
    counts[alphabet[0]] = 0;
    if(I.size() == 0){
        for(int64_t i = 1; i < alphabet.size(); i++)
            counts[alphabet[i]] = 0;
        return;
    }
    
    for(int64_t i = 1; i < alphabet.size(); i++){
        int64_t count_prev = bwt.rank(I.right+1, alphabet[i-1]) - bwt.rank(I.left, alphabet[i-1]) ;
        counts[alphabet[i]] = counts[alphabet[i-1]] + count_prev;
    }
    
}

template<class t_bitvector>
bool is_right_maximal(const BD_BWT_index<t_bitvector>& index, Interval_pair I){
    
    // An interval is right-maximal iff it has more than one possible right extension
    vector<uint8_t> symbols = get_interval_symbols(index.reverse_bwt, I.reverse);
    return (symbols.size() >= 2 || (I.reverse.size() >= 2 && symbols[0] == '$')); // TODO: don't hardcode the dollar
}

template<class t_bitvector>
bool is_left_maximal(const BD_BWT_index<t_bitvector>& index, Interval_pair I){ // UNTESTED
    
    // An interval is left-maximal iff it has more than one possible left extension
    vector<uint8_t> symbols = get_interval_symbols(index.forward_bwt, I.forward);
    return (symbols.size() >= 2 || (I.forward.size() >= 2 && symbols[0] == '$')); // TODO: don't hardcode the dollar
}

// Returns the alphabet sorted order
vector<uint8_t> get_alphabet(string& s){
    
    vector<bool> found(256,false);
    for(uint8_t c : s) found[c] = true;
    vector<uint8_t> alphabet;
    for(int i = 0; i < 256; i++){
        if(found[i]) alphabet.push_back((char)i);
    }
    
    return alphabet;
}

template<class t_bitvector>
void initialize_index_external_memory_psascan(BD_BWT_index<t_bitvector>& index, string& input, string workspace_path, string psascan_path){
    // currently should have 8n memory (the input string)
    
    index.alphabet = get_alphabet(input);
    cerr << getTimeString() << " Building the BWT of the forward string" << endl;
    string forward_bwt_path = workspace_path + "/forward.bwt";
    bwt_pSAscan(input, workspace_path, psascan_path, forward_bwt_path, input.size() * 3 / 2);
    write_to_disk(input, workspace_path + "/inputstring_temp.txt"); input = ""; input.shrink_to_fit(); // Put input to disk
    
    // 0 memory
    cerr << getTimeString() << " Building a wavelet tree out of the forward BWT" << endl;
    construct(index.forward_bwt, forward_bwt_path, 1); // [wt construction] memory
    // [wt] memory
    input = read_from_disk(workspace_path + "/inputstring_temp.txt"); // Read input from disk
    // [wt] + 8n memory
    // Build backward bwt
    std::reverse(input.begin(), input.end());
    
    cerr << getTimeString() << " Building the BWT of the backward string" << endl;
    string backward_bwt_path = workspace_path + "/backward.bwt";
    bwt_pSAscan(input, workspace_path, psascan_path, backward_bwt_path, input.size() * 3 / 2);
    // [wt] + 8n memory
    std::reverse(input.begin(), input.end());
    
    // Build wavelet tree out of backward bwt
    input = ""; input.shrink_to_fit(); // Free input from memory
    // [wt] memory]
    cerr << getTimeString() << " Building a wavelet tree out of the backward BWT" << endl;
    construct(index.reverse_bwt, backward_bwt_path, 1);  // [wt] + [wt construction] memory
    // 2*[wt] memory
    // Don't load input back anymore
    
    // Compute cumulative character counts
    index.cumulative_char_count.resize(256);
    count_smaller_chars(index.forward_bwt,index.alphabet,index.cumulative_char_count,Interval(0,index.forward_bwt.size()-1));
    
    assert(index.forward_bwt.sigma == index.alphabet.size()); // Make sure no newlines or nulls are in the wavelet tree
}

template<class t_bitvector>
void initialize_index_external_memory_ropebwt(BD_BWT_index<t_bitvector>& index, string input_fasta_path, string workspace_path, int64_t extra_memory_bytes){

    string forward_sorted_rotated_path = workspace_path + "/forward_sorted_rotated.fna";
    
    cerr << getTimeString() << " Preprocessing for ropebwt" << endl;
    ropebwt_preprocess(input_fasta_path, forward_sorted_rotated_path, false); // Sort and rotate
    
    cerr << getTimeString() << " Building the BWT of the forward string (ropebwt)" << endl;
    string forward_bwt_path = workspace_path + "/forward.bwt";
    bwt_ropebwt(forward_sorted_rotated_path, forward_bwt_path, workspace_path, extra_memory_bytes);
    cerr << getTimeString() << " Building a wavelet tree out of the forward BWT" << endl;
    construct(index.forward_bwt, forward_bwt_path, 1);
    
    cerr << getTimeString() << " Preprocessing for reverse ropebwt" << endl;
    string reverse_path = workspace_path + "/reverse_fasta.fna"; 
    reverse_fasta(input_fasta_path, reverse_path);

    cerr << getTimeString() << " Building the BWT of the backward string (ropebwt)" << endl;
    
    // NOTE: The resulting BWT is not correct for the positions corresponding to suffixes
    // that start with a dollar. This means we can't invert the reverse BWT, but we can
    // backward search it, as long as we don't search into dollars.
    string backward_bwt_path = workspace_path + "/backward.bwt";
    bwt_ropebwt(reverse_path, backward_bwt_path, workspace_path, extra_memory_bytes);
    
    // Build wavelet tree out of backward bwt
    cerr << getTimeString() << " Building a wavelet tree out of the backward BWT" << endl;
    construct(index.reverse_bwt, backward_bwt_path, 1);  // [wt] + [wt construction] memory

    index.alphabet.clear();
    for(char c : get_interval_symbols(index.forward_bwt, Interval(0,index.forward_bwt.size()-1)))
        index.alphabet.push_back(c);
    sort(index.alphabet.begin(), index.alphabet.end());
    
    // Compute cumulative character counts
    index.cumulative_char_count.resize(256);
    count_smaller_chars(index.forward_bwt,index.alphabet,index.cumulative_char_count,Interval(0,index.forward_bwt.size()-1));
    
}

template<class t_bitvector>
void initialize_index_internal_memory(BD_BWT_index<t_bitvector>& index, string& input){
    index.alphabet = get_alphabet(input);
    
    char* forward = (char*) malloc(input.size() + 1);
    char* backward = (char*) malloc(input.size() + 1);
    for(int64_t i = 0; i < input.size(); i++){
        forward[i] = input[i];
        backward[input.size()-1-i] = input[i];
    }
    forward[input.size()] = 0;
    backward[input.size()] = 0;

    char* forward_bwt = bwt_dbwt(forward);
    char* backward_bwt = bwt_dbwt(backward);
    free(forward);
    free(backward);
    
    // Build wavelet trees
    construct_im(index.forward_bwt, string(forward_bwt).c_str(), 1); // Not casting through std::string segfaults
    construct_im(index.reverse_bwt, string(backward_bwt).c_str(), 1); // Not casting through std::string segfaults
    free(forward_bwt);
    free(backward_bwt);
    
    // Compute cumulative character counts
    index.cumulative_char_count.resize(256);
    count_smaller_chars(index.forward_bwt,index.alphabet,index.cumulative_char_count,Interval(0,index.forward_bwt.size()-1));

}

template<class t_bitvector>
bool allDistinct(const wt_huff<t_bitvector>& bwt, Interval I){
    return get_interval_symbols<t_bitvector>(bwt,I).size() == I.size();
}

template<class t_bitvector>
bool interval_is_supermaximal(const BD_BWT_index<t_bitvector>& index, Interval_pair I){
    return allDistinct(index.forward_bwt,I.forward) && allDistinct(index.reverse_bwt,I.reverse);
}

