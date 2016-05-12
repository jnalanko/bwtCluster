#ifndef SDSL_ITERATE_HH
#define SDSL_ITERATE_HH

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include <vector>
#include <utility>
#include <deque>
#include <string>
#include "Interval.hh"

// Bidirectional BWT index
template<class t_bitvector>
class BD_BWT_index{
private:
    
public:
    
    sdsl::wt_huff<t_bitvector> forward_bwt;
    sdsl::wt_huff<t_bitvector> reverse_bwt;
    
    // cumulative_char_count[i] = total number of characters in the string with ascii value strictly less than i
    std::vector<int64_t> cumulative_char_count;
    
    std::vector<uint8_t> alphabet; // TODO: reverse complement alphabet?

};


/**
 * @brief    Function initialize_index_external_memory_psascan
 * Initializes the given index using external memory and psascan. NOTE: Will clear given the input string
 *
 * @param    index          The index to be initialized
 * @param    input          The input data string (WILL BE CLEARED)
 * @param    workspace_path The path of the working space directory
*/
template<class t_bitvector>
void initialize_index_external_memory_psascan(BD_BWT_index<t_bitvector>& index, std::string& input, std::string workspace_path, std::string psascan_path);


template<class t_bitvector>
void initialize_index_external_memory_ropebwt(BD_BWT_index<t_bitvector>& index, std::string input_fasta_path, 
                                              std::string workspace_path, int64_t extra_memory_megabytes);


/**
 * @brief    Function initialize_index_internal_memory
 *
 * Initializes the given index.
 *
 * @param    index    The index to be initialized
 * @param    input    The input data string
*/
template<class t_bitvector>
void initialize_index_internal_memory(BD_BWT_index<t_bitvector>& index, std::string& input);

template<class t_bitvector>
void count_smaller_chars(const sdsl::wt_huff<t_bitvector>& bwt, std::vector<uint8_t>& alphabet, std::vector<int64_t>& counts, Interval I);

// Returns an interval of size zero extension not possible or if the given interval has size 0
template<class t_bitvector>
Interval_pair left_extend(const BD_BWT_index<t_bitvector>& index, Interval_pair intervals, char c);

// Returns an interval of size zero extension not possible or if the given interval has size 0
template<class t_bitvector>
Interval_pair left_extend(const BD_BWT_index<t_bitvector>& index, Interval_pair intervals, char c, int64_t cumul_rank_c);

// Returns an interval of size zero extension not possible or if the given interval has size 0
template<class t_bitvector>
Interval_pair right_extend(const BD_BWT_index<t_bitvector>& index, Interval_pair intervals, char c);

// Returns an interval of size zero extension not possible or if the given interval has size 0
template<class t_bitvector>
Interval_pair right_extend(const BD_BWT_index<t_bitvector>& index, Interval_pair intervals, char c, int64_t cumul_rank_c);

template<class t_bitvector>
int64_t backward_step(const BD_BWT_index<t_bitvector>& index, int64_t pos);

template<class t_bitvector>
bool interval_is_supermaximal(const BD_BWT_index<t_bitvector>& index, std::pair<Interval,Interval> I);


// Returns the symbols in the interval as a vector in arbitary order
template<class t_bitvector>
std::vector<uint8_t> get_interval_symbols(const sdsl::wt_huff<t_bitvector>& wt, Interval I);

// Another version of get_interval_symbols when the allocation of a new vector is too slow.
// Counts the number of distinct symbols into nExtensions, and puts the symbols into the indices
// [0,nExtensions[ of symbols. Also stores ranks of the symbols at the endpoints of the interval I
// to ranks_i and ranks_j. Important: All the parameter vectors must have length at least equal to the size
// of the alphabet of the given wavelet tree. Symbols is not sorted
template<class t_bitvector>
void get_interval_symbols(const sdsl::wt_huff<t_bitvector>& wt, Interval I, sdsl::int_vector_size_type& nExtensions, std::vector<uint8_t>& symbols,
  std::vector<uint64_t>& ranks_i, std::vector<uint64_t>& ranks_j);

template<class t_bitvector>
int64_t compute_cumulative_char_rank_in_interval(const sdsl::wt_huff<t_bitvector>& wt, std::vector<char> alphabet, char c, Interval I);

#include "BD_BWT_index_impl.hh"

#endif
