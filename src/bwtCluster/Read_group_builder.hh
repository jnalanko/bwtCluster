#ifndef READ_GROUP_BUILDER_HH
#define READ_GROUP_BUILDER_HH

#include "UnionFind.hh"

template <class t_bitvector>
class Read_group_builder{
private:
    BD_BWT_index<t_bitvector>* index;
    int64_t nReads;
    int64_t nThreads;
    std::vector<BWT_inversion_iterator<t_bitvector> > iterators;
    std::mutex union_find_mutex;
    std::mutex interval_reps_mutex;
    std::string paired_end_file;
    
    void write_groups(UnionFind& UF, std::ostream &output);
    sdsl::bit_vector compute_k_submaximal_intervals(int64_t k);
    UnionFind parallel_non_RC_union_find(int64_t maxGroupSize, sdsl::bit_vector& intervals);
    void do_unions(BWT_inversion_iterator<t_bitvector>iter, sdsl::int_vector<0>& interval_reps, 
                   UnionFind& UF, sdsl::bit_vector& intervals, sdsl::rank_support_v<1>& intervals_rs, 
                   int64_t maxGroupSize);
    void compute_marked_position_read_ids(BWT_inversion_iterator<t_bitvector> iter, sdsl::bit_vector& marks,
            sdsl::rank_support_v<1>& marks_rs, sdsl::int_vector<0>& read_ids, std::mutex& mutex); 
public:
    
     /**
     * @brief    Constructor Read_group_builder
     * 
     * Constructs a Read group builder that uses as many threads as there are iterators in the parameter vector
     * 
     * @param    index        The bidirectional index
     * @param    iterators    BWT inversion iterators for the parallel inversion of the BWT
     */
    Read_group_builder<t_bitvector>(BD_BWT_index<t_bitvector>* index, std::vector<BWT_inversion_iterator<t_bitvector> > iterators, std::string paired_end_file) 
        : index(index), iterators(iterators), paired_end_file(paired_end_file){
        // Number of reads: One dollar for each read + one dollar at the end of data
        nReads = index->forward_bwt.rank(index->forward_bwt.size(),'$') - 1;
        nThreads = iterators.size();
    }
    
    /**
     * @brief Function groupReads
     *
     * Groups two reads together if they share a common k-mer. The grouping is written to the given
     * output stream such that each line contains space delimited read IDs of one group
     *
     * @param    k                The k-mer length
     * @param    maxGroupSize    The maximum size of a group
     * @param    output            A stream to write the output into
     */
    void groupReads(int64_t k, int64_t maxGroupSize, std::ostream &output);
    
    /**
     * @brief Function groupReads
     *
     * Groups two reads together if they share a common k-mer. Implemented using k-submaximal repeats to save space. 
     * The grouping is written to the given output stream such that each line contains space delimited read IDs of one group
     *
     * @param    k                The k-mer length
     * @param    maxGroupSize    The maximum size of a group
     * @param    output            A stream to write the output into
     */
    void groupReadsSubmaximal(int64_t k, int64_t maxGroupSize, std::ostream &output);
    
    void groupReadsRC(int64_t k, int64_t maxGroupSize, std::string workspace_path, std::ostream &output);
    
};

#include "Read_group_builder_impl.hh"

#endif