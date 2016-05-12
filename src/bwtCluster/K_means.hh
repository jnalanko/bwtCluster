#ifndef K_MEANS_H_
#define K_MEANS_H_

#include <iostream>
#include <vector>
#include <string>

// Metacluster's implementation of k-means
class K_means_metacluster {
private:
    std::string precluster_file;
    std::string out_dir;
    int64_t kmer_k;
public:

    K_means_metacluster(std::string precluster_file, std::string out_dir, int64_t kmer_k) :
        precluster_file(precluster_file), out_dir(out_dir), kmer_k(kmer_k) {}
    virtual int64_t run(std::string results_file, std::string workspace); // Return the number of clusters

};

// Our implementation of k-means
// NOTE: Assumes the alphabet {A,C,G,T}
class K_means {

private:

    std::string precluster_file; // Filename of the file with the preclusters
    std::vector<std::vector<float> > distributions; // kmer distributions of preclusters
    std::vector<std::vector<int32_t> > rankVectors; // rank vectors of kmer distributions of preclusters

    // Compute the distribution of k_short-mers that are covered by a k_long-mer that occurs at least tau times
    std::vector<float> processGroup(std::vector<std::string>& precluster, int64_t k_long, int64_t k_short, bool rc, int64_t tau);

    // Compute the lexicographical rank of kmer among all k-length string with alphabet {A,C,G,T}
    int64_t lexRank(std::string kmer);
    
    // Computes the rank vector of the given distribution. The rank of a position
    // is the number of positions that have a value less than that position. Ties
    // are broken arbitrarily
    std::vector<int32_t> computeRankVector(const std::vector<float>& distribution);
    
    // Computes the closest centroid to the given rank vector point using the Spearman footrule distance
    int64_t getBestCentroid(const std::vector<int32_t>& point, const std::vector<std::vector<int32_t> >& centroids);
    
    // Finds the complements character of c, i.e. A -> T, C -> G, G -> C, T -> A
    char complement(char c);
    
    // Find the reverse-complement string of s
    std::string getRC(std::string s);
    
    // Writes clusters out with one line per cluster. Filters preclusters that
    // are at least two standard deviations from the mean distance from the center of the cluster
    void write_clusters_without_outliers(std::vector<std::vector<float> >& centroids, 
                                         std::vector<std::vector<int64_t> >& assignments,
                                         std::vector<std::vector<int32_t> >& centroidRanks,
                                         std::string results_directory, std::string workspace_dir);

    double interClusterDistance(std::vector<std::vector<int64_t> >& clusters, int64_t c1, int64_t c2);

    std::vector<std::vector<int64_t> > mergeClusters(double threshold, std::vector<std::vector<int64_t> > & assignments);
    void writeClusters(std::vector<std::vector<int64_t> >& finalClusters, std::string results_file, std::string workspace_dir);

    
public:


    /** @brief Constructor
     *
     * Loads the preclusters into memory and computes the kmer distributions
     * of the frequent repeats of the reads. The preclusters are given in a fasta file,
     * where the preclusters are consecutive reads in the file such that the header of the
     * first read of each precluster is ">precluster_start"
     *
     * @param    precluster_file  The fasta file which contains the preclusters.
     * @param    rc        True = consider reverse complements
     * @param    kmer_k    The k-mer length
     *
     */
    K_means(std::string precluster_file, bool rc, int64_t kmer_k);

    /**
     * @brief Function clusterGroups
     *
     * Clusters groups based on the Spearman footrule distance of the kmer compositions
     * and the k-means algorithm
     *
     * @param    kmeans_k    The parameter k of k-means
     * @param    nRounds     Number of rounds to run the k-means algorithm
     * @param    results_file The filename of the output
     * @param    workspace    Working space directory
     * 
     *
     */
    void run(int64_t kmeans_k, int64_t nRounds, std::string results_file, std::string workspace);

};
#endif /* K_MEANS_H_ */
