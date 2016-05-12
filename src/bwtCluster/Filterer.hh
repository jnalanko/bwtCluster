#ifndef FILTERER_HH
#define FILTERER_HH

#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "Iterators.hh"

template <class t_bitvector>
class Filterer{
private:
    BD_BWT_index<t_bitvector>* index;
    int64_t nThreads;
    std::mutex output_mutex;
    std::vector<BWT_inversion_iterator<t_bitvector> > iterators;
    std::string workspace_path;
    std::string paired_end_file;
    
    void output_marked_reads(BWT_inversion_iterator<t_bitvector> iter, std::ostream &goodReads, std::ostream &badReads, 
                         std::vector<int64_t>& goods_ids, const sdsl::bit_vector& marks);
    
public:
     /**
     * @brief    Constructor Filterer
     *
     * Constructs a Filterer that uses as many threads as there are iterators in the parameter vector
     * 
     * @param index           The bidirectional index
     * @param iterators       BWT inversion iterators for the parallel inversion of the BWT
     * @param workspace_path  A path to a working space directory, without a trailing slash
     * @param paired_end_file The path to the file that contains the ranks of the paired reads, one pair per line
     */    
    Filterer(BD_BWT_index<t_bitvector>* index, std::vector<BWT_inversion_iterator<t_bitvector> > iterators, std::string workspace_path, std::string paired_end_file);

     /**
     * @brief    Function filterReadsRC
     *
     * Filter reads with no k-mer with at least tau occurrences. TODO: update
     *
     * @param    good_reads    A stream to write the good reads into in FASTA format
     * @param    bad_reads     A stream to write the bad reads into in FASTA format
     * @param    good_ids      A stream to write ids of the kept reads, one on each line
     * @param    bad_ids       A stream to write ids of the not kept reads, one on each line
     * 
     * @return Ranks of reads that were kept
     */
    vector<int64_t> filterReads(int k, int tau, std::string good_out, std::string bad_out);
    
     /**
     * @brief    Function filterReadsRC
     *
     * Filter reads with no k-mer with at least tau occurrences (RC). TODO: update
     *
     * @param    good_reads    A stream to write the good reads into in FASTA format
     * @param    bad_reads     A stream to write the bad reads into in FASTA format
     * @param    good_ids      A stream to write ids of the kept reads, one on each line
     * @param    bad_ids       A stream to write ids of the not kept reads, one on each line
     * 
     * @return Ranks of reads that were kept
     */
     
    vector<int64_t> filterReadsRC(int k, int tau, std::string good_out, std::string bad_out);

};

#include "Filterer_impl.hh"

#endif