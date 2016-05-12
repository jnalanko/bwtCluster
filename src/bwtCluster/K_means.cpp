#include "K_means.hh"
#include "Parser.hh"
#include "tools.hh"
#include "UnionFind.hh"
#include <iostream>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <set>
#include <fstream>
#include <cassert>
#include <thread>
#include <unordered_map>

using namespace std;

static int64_t power(int64_t base, int64_t power) {
    int64_t result = 1;
    for(int64_t i = 0; i < power; i++){
        result *= base;
    }
    return result;
}

int64_t K_means::lexRank(string kmer){
    // Compute the rank by interpreting the kmer as a number
    // written in base 4 with A=0, C=1, G=2, T=3
    int64_t rank = 0;
    int64_t power_of_4 = 1;
    for(int64_t i = kmer.size() - 1; i >= 0; i--){
        char c = kmer[i];
        if(c == 'a' || c == 'A') rank += power_of_4 * 0; // TODO: do we need the lower case?
        if(c == 'c' || c == 'C') rank += power_of_4 * 1;
        if(c == 'g' || c == 'G') rank += power_of_4 * 2;
        if(c == 't' || c == 'T') rank += power_of_4 * 3;
        power_of_4 *= 4;
    }
    return rank + 1;
}

static string unrank(int64_t x, int64_t k){
    // Base-4 representation of x with 0 = A, 1 = C, 2 = G, 3 = C
    x--; // Switch to zero indexing
    string ans = "";
    int64_t b = power(4,k-1);
    while(b > 0){
        int64_t digit = x / b;
        switch(digit){
            case 0: ans += 'A'; break;
            case 1: ans += 'C'; break;
            case 2: ans += 'G'; break;
            case 3: ans += 'T'; break;
        }
        x -= (x/b)*b;
        b /= 4;
    }
    return ans;
}

vector<int32_t> K_means::computeRankVector(const vector<float>& distribution){
    vector<pair<double, int64_t> > v;
    for(int64_t i = 0; i < distribution.size(); i++){
        v.push_back(make_pair(distribution[i],i));
    }
    sort(v.begin(), v.end()); // Sort by the first element of the pair, i.e. the distribution value
    vector<int32_t> result(distribution.size());
    for(int64_t i = 0; i < distribution.size(); i++){
        int64_t pos = v[i].second;
        int64_t rank = i;
        result[pos] = rank;
    }
    return result;
}
    
char K_means::complement(char c){
    if(c == 'A') return 'T';
    if(c == 'C') return 'G';
    if(c == 'G') return 'C';
    if(c == 'T') return 'A';
    cerr << "Invalid character: '" << c << "'" << endl; 
    assert(false);
}

string K_means::getRC(string s){
    string ans = "";
    for(int i = s.size() - 1; i >= 0; i--){
        ans += complement(s[i]);
    }
    return ans;
}

vector<float> K_means::processGroup(vector<string>& precluster, int64_t k_long, int64_t k_short, bool rc, int64_t tau){

    unordered_map<string,int64_t> kmers;
    
    // Count the k_long-mers in all reads
    for(string& read : precluster){
        if(read.size() < k_long) cout << read << endl;
        assert(read.size() >= k_long);
        for(int64_t pos = 0; pos <= read.size() - k_long; pos++){
            kmers[read.substr(pos,k_long)]++;
        }
    }
    
    
    vector<float> dist(power(4,k_short),0);
    int64_t nKmers = 0;
    
    // Compute the k_short-mer distribution of k_short_mers inside the k_long-mers that have frequency
    // at least tau.
    for(auto keyvalue : kmers){
        string kmer = keyvalue.first;
        int64_t count = keyvalue.second;
        if(count >= tau){
            for(int64_t i = 0; i <= k_long - k_short; i++){
                dist[lexRank(kmer.substr(i,k_short))]++;
                nKmers++;
            }
        }
    }
    
    if(rc){
        // Take the lexicographically smallest of a reverse-complement pair as a
        // representative for that kmer
        vector<float> dist_rc;
        for(int64_t i = 1; i <= dist.size(); i++){
            string kmer = unrank(i, k_short);
            string kmer_rc = getRC(kmer);
            if(lexRank(kmer) <= lexRank(kmer_rc)){
                if(lexRank(kmer) < lexRank(kmer_rc)) 
                    dist_rc.push_back(dist[lexRank(kmer)] + dist[lexRank(kmer_rc)]);
                else
                    dist_rc.push_back(dist[lexRank(kmer)]);
            }
        }
        dist = dist_rc;
    }
        
    for(int64_t i = 0; i < dist.size(); i++){
        dist[i] /= nKmers; // Normalize
    }
    
    return dist;
   
}


void K_means::worker_thread(ParallelQueue<std::pair<std::vector<string>, int64_t> >& Q, int64_t kmer_k, bool rc){
    while(true){
        pair<vector<string>,int64_t> work_frame = Q.pop();
        vector<string> precluster = work_frame.first;
        int64_t id = work_frame.second;
        if(id == -1){
            // No more work signal received
            Q.push(work_frame);
            break;
        }
        vector<float> distribution = processGroup(precluster, 16, kmer_k, rc, 2); // TODO: 16 and 2 into the config file
        
        // Store the results
        distributions[id] = distribution;
        rankVectors[id] = computeRankVector(distribution);
    }    
}


K_means::K_means(string precluster_file, bool rc, int64_t kmer_k, int nThreads){
    if(rc)
        cerr << getTimeString() << " Computing " << kmer_k << "-mer distributions for k-means (RC)" << endl;
    else
        cerr << getTimeString() << " Computing " << kmer_k << "-mer distributions for k-means" << endl;

    this->precluster_file = precluster_file;

    
    Parser parser;
    
    // Count the number of preclusters
    ifstream infile_counting(precluster_file);
    int64_t number_of_preclusters = parser.count_clusters(infile_counting);
    cerr << getTimeString() << " " << number_of_preclusters << " preclusters in k-means" << endl;
    infile_counting.close();
        
    // Allocate space for the distributions and rank vectors
    distributions.resize(number_of_preclusters);
    rankVectors.resize(number_of_preclusters);
    distributions.shrink_to_fit();
    rankVectors.shrink_to_fit();
    
    // A parallel queue to preprocess the preclusters
    // Each element in the queue the set of reads in one precluster, and
    // the id of the precluster.
    // The id -1 in the queue signals the worker threads to quit
    ParallelQueue<std::pair<std::vector<string>, int64_t> > Q;
        
    // Create the threads
    vector<thread> threads(nThreads);
    for(int i = 0; i < nThreads; i++){
        threads[i] = std::thread(&K_means::worker_thread, this, ref(Q), kmer_k, rc);
    }
    
    // Push work to the queue
    ifstream infile(precluster_file);
    int64_t precluster_id = 0;
    while(true){
        vector<string> precluster = parser.next_cluster(infile);
        if(precluster.size() != 0)
            Q.push({precluster,precluster_id});
        else{
            Q.push({precluster,-1});
            break;
        }
        precluster_id++;
    }
    
    assert(precluster_id == number_of_preclusters);
    
    // Wait for the threads to finish
    for(int i = 0; i < nThreads; i++) threads[i].join();

}


static int32_t rankCorrelation(const vector<int32_t> x1, const vector<int32_t> x2){
    int64_t result = 0;
    for(int64_t i = 0; i < x1.size(); i++){
        result += abs(x1[i] - x2[i]);
    }
    return result;
}

double K_means::interClusterDistance(vector<vector<int64_t> > & clusters, int64_t c1, int64_t c2){
    double distance_sum = 0;
    int64_t nPairs = 0;
    for(int64_t x1 : clusters[c1]){
        for(int64_t x2 : clusters[c2]){
            if(x1 == x2) continue; // Happens only if c1 == c2
            distance_sum += rankCorrelation(rankVectors[x1], rankVectors[x2]);
            nPairs++;
        }
    }
    if(nPairs == 0) return 0;
    return distance_sum / nPairs;
}

vector<vector<int64_t> >  K_means::mergeClusters(double threshold, vector<vector<int64_t> > & assignments){
    // Remove empty clusters
    vector<vector<int64_t> > clusters;
    for(auto v : assignments)
        if(v.size() != 0)
            clusters.push_back(v);
    
    vector<double> intraDistances;
    for(int64_t i = 0; i < clusters.size(); i++){
        intraDistances.push_back(interClusterDistance(clusters,i,i));
    }

    UnionFind UF;
    UF.init(clusters.size());
    int64_t nJoins = 0;
    for(int64_t i = 0; i < clusters.size(); i++){
        for(int64_t j = i + 1; j < clusters.size(); j++){
            if(clusters.size() - nJoins <= 4) break;
           // cerr << interClusterDistance(clusters,i,j) << " " << (intraDistances[i] + intraDistances[j])/2.0 << endl;
            if(threshold * interClusterDistance(clusters,i,j) <= (intraDistances[i] + intraDistances[j])/2.0){
                if(UF.find(i) != UF.find(j)){
                    UF.doUnion(i,j);
                    nJoins++;
                }
            }
        }
    }
    cerr << getTimeString() << " Final number of clusters: " << clusters.size() - nJoins << endl;
    
    vector<vector<int64_t> > finalClusters(clusters.size());
    for(int64_t i = 0; i < clusters.size(); i++){
        int64_t final_id = UF.find(i);
        for(int64_t x : clusters[i])
            finalClusters[final_id].push_back(x);
    }
    
    // Remove empty clusters
    int64_t j = 0;
    for(int64_t i = 0; i < finalClusters.size(); i++){
        if(finalClusters[i].size() != 0){
            if(i != j) // Avoid assignemtn to itself
                finalClusters[j] = finalClusters[i];
            j++;
        }
    }
    
    finalClusters.resize(j);
    return finalClusters;
}

int64_t K_means::getBestCentroid(const vector<int32_t>& point, const vector<vector<int32_t> >& centroids){
    int64_t bestCentroid = -1;
    int64_t bestCorrelation = 999999999999999LL;
    for(int64_t i = 0; i < centroids.size(); i++){
        int64_t correlation = rankCorrelation(point, centroids[i]);
        {
            if(correlation < bestCorrelation){
                bestCorrelation = correlation;
                bestCentroid = i;
            }
        }
    }
    return bestCentroid;
}

void K_means::writeClusters(vector<vector<int64_t> >& finalClusters, string results_file, string workspace_dir){

    Parser parser;
    ifstream infile(precluster_file);
    vector<int64_t> precluster_to_finalcluster(distributions.size());
    
    int64_t k = finalClusters.size();

    // Open the output streams
    vector<ofstream*> outputs;
    for(int64_t i = 0; i < k; i++){
        ofstream* out = new ofstream(workspace_dir + "/cluster_" + to_string(i));
        if(!out->good()){
            cerr << "Error opening file " << workspace_dir + "/cluster_" + to_string(i) << endl;
            exit(1);
        }
        outputs.push_back(out);
    }
    
    for(int64_t center = 0; center < k; center++){

        // Add all reads of the precluster to the cluster set
        for(int64_t precluster : finalClusters[center]){
            precluster_to_finalcluster[precluster] = center;
        }
    }
    
    // Write the final clusters into files
    for(int64_t i = 0; i < distributions.size(); i++){
        vector<string> precluster = parser.next_cluster(infile);
        int64_t center = precluster_to_finalcluster[i];
        for(auto& read : precluster){
            *(outputs[center]) << ">\n" << read << "\n";
        }
    }

    // Clean up
    for(int i = 0; i < k; i++){
        outputs[i]->flush();
        delete outputs[i];
    }
            
    // Concatenate all clusters into one file
    ofstream results(results_file);
    for(int64_t i = 0; i < k; i++){
        ifstream in(workspace_dir + "/cluster_" + to_string(i));
        int64_t count = 0;
        while(true){
            pair<string,string> x = parser.next_read(in);
            string read = x.first;
            if(x.first == "" && x.second == "") break;
            
            if(count == 0) results << ">cluster_start\n";
            else results << ">\n";
            
            results << read << "\n";
            count++;
        }
        in.close();
        
        // Clear the file
        ofstream X(workspace_dir + "/cluster_" + to_string(i));
        X.close();
    }
    
    results.close(); // flush

}


void K_means::run(int64_t kmeans_k, int64_t nRounds, string results_file, string workspace){
    // Assignment phase: iterate through all reads, compute the distance to all centroids, assign to closest
    // Update phase: for all reads in a cluster, take the arithmetic mean.
    
    assert(distributions.size() >= kmeans_k);
    cerr << getTimeString() << " Running k-means" << endl;

    if(distributions.size() == 0){
        throw runtime_error("Error: K-means called with no input loaded");
    }

    // Initialize with random centroids
    srand(time(NULL));
    vector<vector<float> > centroids;
    int64_t distLength = distributions[0].size(); // Distribution vector length
    vector<int64_t> random_indices;
    for(int64_t i = 0; i < kmeans_k; i++){
        int64_t pos = rand() % distributions.size();
        while(find(random_indices.begin(), random_indices.end(), pos) != random_indices.end())
            pos = rand() % distributions.size();
        random_indices.push_back(pos);
        vector<float> d = distributions[pos];
        centroids.push_back(d);
    }

    int64_t roundCount = 1;
    vector<vector<int64_t> > prev_assignments;
    while(roundCount <= nRounds){
        // vector sums: For each centroid, the sum of all the distributions assigned to that centroid
        vector<vector<float> > sums(centroids.size());

        // vector assignment: A list of precluster ids for each centroid
        vector<vector<int64_t> > assignments(centroids.size());

        // Initialize and compute the rank vectors for all centroids
        vector<vector<int32_t> > centroidRanks(centroids.size());
        
        #pragma omp parallel for
        for(int i = 0; i < centroids.size(); i++){
            centroidRanks[i] = computeRankVector(centroids[i]);
            sums[i] = vector<float>(distLength,0);
            assignments[i] = vector<int64_t>();
        }

        // Assign groups and accumulate the sum of the distributions for all centroids
        vector<int64_t> nAssignmentsForCentroid(centroids.size(),0);
        for(int64_t i = 0; i < distributions.size(); i++){
            int64_t bestCentroid = getBestCentroid(rankVectors[i], centroidRanks);
            assignments[bestCentroid].push_back(i);
            for(int j = 0; j < distLength; j++){
                sums[bestCentroid][j] += distributions[i][j];
            }
            nAssignmentsForCentroid[bestCentroid] += 1;
        }
        

        // Divide the accumulated sums by the number of groups for each centroid
        #pragma omp parallel for
        for(int64_t centroid = 0; centroid < centroids.size(); centroid++){
            if(nAssignmentsForCentroid[centroid] == 0) continue;
            for(int j = 0; j < distLength; j++){
                sums[centroid][j] /= nAssignmentsForCentroid[centroid];
            }
        }

        // If finished, print the output
        for(auto& v : assignments)
            sort(v.begin(), v.end());
        
        if(roundCount == nRounds || assignments == prev_assignments){
            //vector<vector<int64_t> > finalClusters = mergeClusters(1.00,assignments); // 0.79 from Metacluster 3.0 paper
            writeClusters(assignments,results_file,workspace);
            cerr << getTimeString() << " K-means finished after " << roundCount << " rounds" << endl;
            return;
        }
        
        roundCount += 1;
        centroids = sums;
        prev_assignments = assignments;
    }

}