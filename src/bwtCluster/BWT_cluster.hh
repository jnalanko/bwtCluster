#ifndef BWT_CLUSTER
#define BWT_CLUSTER

#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "Iterators.hh"
#include "bwt.hh"
#include "Filterer.hh"
#include "Algorithms.hh"

class UnionFind;

template <class t_bitvector>
class BWT_cluster{
    
private:
        
    // Construct a bwt of the string using psascan
    void buildBWT_psascan(std::string& input, std::string psascan_path);
    
    // Construct a BWT of the string using ropebwt
    void buildBWT_ropebwt(std::string fasta_filename);
    
    // Write the given space and a description to the memory log along with the current time stamp
    void memlog(int64_t space, std::string label);
    
    // Finds the bwt position of the read by backwards searching the read
    // If the read is not unique, crashes by failing an assertion.
    int64_t locate_read_start(const string& read);
    
    // Get parallel partition assuming the reads are lexicographically sorted in the concatenation.
    std::vector<BWT_inversion_iterator<t_bitvector> > get_iterators_sorted();
    
    BD_BWT_index<t_bitvector> index; // The bidirectional BWT index of the concatenation of the reads
    int64_t dataLength; // Number of characters in concatenation string
    int64_t nReads; // Number of reads in data
    int64_t space_bits; // Number of bits of space used at the moment in theory
    int64_t start_time; // The time unix epoch in milliseconds when the object was constructed
    int64_t nThreads; // Number of parallel threads to use
    std::string workspace_path; // The path to a directory in disk which will be used as working space
    std::string paired_end_file; // The path to the file that contains the ranks of the paired reads, one pair per line
    ofstream memory_log; // The stream to write the memory logging messages
    
public:
    
    /**
     * @brief BWT_cluster constructor
     *
     * @param workspace_path  The directory the object will use for working space
     * @param paired_end_file The path to the file that contains the ranks of the paired reads, one pair per line
     * @param nThreads        Number of parallel threads to use
     */
    BWT_cluster<t_bitvector>(std::string workpace_path, std::string paired_end_file, int64_t nThreads);
    
    ~BWT_cluster<t_bitvector>();

    /**
     * @brief Load data
     *
     * Load a data in FASTA format into the class. The function builds the BWT of the data and computes
     * a read array, which specifies for each position in the BWT the corresponding read ID in the text.
     *
     * @param    input    The input file path
     */
    void loadData(std::string fasta_filename); // Reads data in FASTA format to the object

    /**
     * @brief Filter reads
     *
     * Filter reads with no k-mer with at least tau occurrences. TODO: update
     *
     * @param    k             The length of the k-mers
     * @param    tau           The occurrence count threshold
     * @param    rc            True if should consider reverse complements of kmers
     * @param    good_reads    A stream to write the good reads into in FASTA format
     * @param    bad_reads     A stream to write the bad reads into in FASTA format
     * @param    good_ids      A stream to write ids of the kept reads, one on each line
     * @param    bad_ids       A stream to write ids of the not kept reads, one on each line
     * 
     * @return Ranks of reads that were kept
     */
    vector<int64_t> filterReads(int k, int tau, bool rc, std::string good_reads, std::string bad_reads);

    /**
     * @brief Build preclusters using right-maximal kmers
     *
     * Groups two reads together if they share a common k-mer. The grouping is written to the given
     * output stream such that each line contains space delimited read IDs of one group
     *
     * @param    k                The k-mer length
     * @param    maxGroupSize     The maximum size of a group
     * @param    rc               True if should consider reverse complements of kmers
     * @param    output           A stream to write the output into
     */
    void groupReads(int64_t k, bool rc, int64_t maxGroupSize, std::ostream &output);
    
    /**
     * @brief Build preclusters using submaximal repeats
     *
     * Groups two reads together if they share a common k-mer. Implemented using k-submaximal repeats to save space. 
     * The grouping is written to the given output stream such that each line contains space delimited read IDs of one group
     *
     * @param    k               The k-mer length
     * @param    maxGroupSize    The maximum size of a group
     * @param    rc              True if should consider reverse complements of kmers. NOTE: not implemented for rc == true
     * @param    output          A stream to write the output into
     */
    void groupReadsSubmaximal(int64_t k, bool rc, int64_t maxGroupSize, std::ostream &output);
    
    /**
     * @brief Filter groups that are too small
     *
     * Discards all groups that have less than a given number of reads. Writes the new groups
     * into the output stream such that each line contains space delimited read IDs of one group
     *
     * @param    groups        The groups in the format outputted by groupReads
     * @param    newGroups     A stream to write the new groups into
     * @param    minSize       All groups that have a size less than this are filtered out
     */
    void filterGroups(std::istream &groups, int64_t minSize, std::ostream &newGroups);

};

#include "BWT_cluster_impl.hh"

#endif


