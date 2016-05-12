#ifndef ALGORITHMS_HH
#define ALGORITHMS_HH

#include <iostream>

/** @brief Function get_right_maximal_kmers
*
* Computes the intervals of the right-maximal k-mers with frequency greater or equal to tau
* for the given index
*
* @param  index    The index
* @param  k        The k-mer length
* @param  tau      The required k-mer frequecy
* @param  fill     If true, mark all positions in the intervals, if false, mark only end points
* @param  nThreads Number of parallel threads to use
*/
template <class t_bitvector>
sdsl::bit_vector get_right_maximal_kmers(BD_BWT_index<t_bitvector>& index, int64_t k, int64_t tau, bool fill, int64_t nThreads);

/** @brief Function get_right_maximal_kmers_rc
*
* Compute the intervals of the right-maximal k-mers with frequency greater or equal to tau
* for the given index, when the reverse complement of a substring is considered the 
* same string as the corresponding forward substring. 
* Returns a bitvectors with all forward suffixes that start with a right-maximal
* kmer marked with ones.
*
* @param  index    The index
* @param  k        The k-mer length
* @param  tau      The required k-mer frequecy
* @param  nThreads Number of parallel threads to use
*/
template <class t_bitvector>
sdsl::bit_vector get_right_maximal_kmers_rc(BD_BWT_index<t_bitvector>& index, int64_t k, int64_t tau, int64_t nThreads);

/** @brief Function find_kmer_rc_pairs
 * 
 * Finds all kmers that appear both in the forward and the RC text. For each such
 * kmer, writes to the given output stream a pair of positions in the forward bwt with an arbitary
 * position from the interval of the forward kmer, and an arbitrary position
 * from the interval of the RC kmer. Returns a bit vector such that position i contains a one
 * if and only if that position was written to disk at some point.
 * 
 * @param index    The index
 * @param k        The k-mer length
 * @param nThreads The number of parallel threads to use
 * @param out      The output stream for the position pairs
 */
template <class t_bitvector>
sdsl::bit_vector find_kmer_rc_pairs(BD_BWT_index<t_bitvector>& index, int64_t k, int64_t nThreads, std::ostream& out);

/** @Brief Function fill_intervals_length_at_least_tau
 * 
 * Fills all intervals with length at least tau with ones.
 * Removes intervals with length less than tau.
 * The input vector marks the start and end position of every
 * interval. It is assumed that the length of every interval is
 * at least 2.
 * 
 * @param v    A bit vector with an even number of ones
 * @param t    The minimum length of an interval to be filled
 * 
 * */
void fill_intervals_length_at_least_tau(sdsl::bit_vector& v, int64_t tau);


#include "Algorithms_impl.hh"

#endif