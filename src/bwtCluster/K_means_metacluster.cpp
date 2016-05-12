#include "K_means.hh"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cassert>
#include <fstream>
#include "MetaCluster.h"
#include "Parser.hh"
#include "tools.hh"

using namespace std;

int64_t power(int64_t base, int64_t power) {
    int64_t result = 1;
    for(int64_t i = 0; i < power; i++){
        result *= base;
    }
    return result;
}

int64_t lexRank(string kmer){
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

string unrank(int64_t x, int64_t k){
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

vector<int64_t> computeRankVector(const vector<int64_t>& distribution){
    vector<pair<int64_t, int64_t> > v;
    for(int64_t i = 0; i < distribution.size(); i++){
        v.push_back(make_pair(distribution[i],i));
    }
    sort(v.begin(), v.end()); // Sort by the first element of the pair, i.e. the distribution value
    vector<int64_t> result(distribution.size());
    for(int64_t i = 0; i < distribution.size(); i++){
        int64_t pos = v[i].second;
        int64_t rank = i;
        result[pos] = rank;
    }
    return result;
}

char complement(char c){
    if(c == 'A') return 'T';
    if(c == 'C') return 'G';
    if(c == 'G') return 'C';
    if(c == 'T') return 'A';
    cerr << "Invalid character: '" << c << "'" << endl; 
    assert(false);
}

string getRC(string s){
    string ans = "";
    for(int i = s.size() - 1; i >= 0; i--){
        ans += complement(s[i]);
    }
    return ans;
}

vector<int64_t> processGroup(vector<string>& precluster, int64_t k_long, int64_t k_short, bool rc, int64_t tau){

    unordered_map<string,int64_t> kmers;
    
    // Count the k_long-mers in all reads
    for(string& read : precluster){
        if(read.size() < k_long) cout << read << endl;
        assert(read.size() >= k_long);
        for(int64_t pos = 0; pos <= read.size() - k_long; pos++){
            kmers[read.substr(pos,k_long)]++;
        }
    }
    
    
    vector<int64_t> dist(power(4,k_short),0);
    int64_t nKmers = 0;
    
    // Count all k_short-mers inside the k_long-mers that have frequency
    // at least tau.
    for(int64_t i = 0; i < precluster.size(); i++){
        string& read = precluster[i];
        int64_t pos_k_short = 0;
        for(int64_t pos_k_long = 0; pos_k_long <= read.size() - k_long; pos_k_long++){
            int64_t k_long_end = pos_k_long + k_long - 1;
            if(kmers[read.substr(pos_k_long,k_long)] >= tau){
                // Add k_short-mers inside the k_long-mer to the distribution
                while(pos_k_short + k_short - 1 <= k_long_end){
                    if(pos_k_short >= pos_k_long){
                        dist[lexRank(read.substr(pos_k_short,k_short))]++;
                        nKmers++;
                    }
                    pos_k_short++;
                }
            }
        }
    }
    
    if(rc){
        for(int64_t i = 1; i <= dist.size(); i++){
            string kmer = unrank(i, k_short);
            string kmer_rc = getRC(kmer);
            if(lexRank(kmer) <= lexRank(kmer_rc)){
                dist[lexRank(kmer)] += dist[lexRank(kmer_rc)];
                dist[lexRank(kmer_rc)] = dist[lexRank(kmer)];
            }
        }
    }
        
    return dist;
   
}
 
int64_t K_means_metacluster::run(string results_file, string workspace){

    Parser parser;
    ifstream infile(precluster_file);
    int64_t nPreclusters = 0;
    int buffer_size = 1024; // Initial buffer, will be doubled if full
    
    int GenoNum = 1; // For some reason always 1
    int** distributions = (int**)malloc(sizeof(int*) * buffer_size); // todo: int64_t. Is now int because Metacluster takes int**
    int** ranks = (int**)malloc(sizeof(int*) * buffer_size); // todo: int64_t
    int** Component = (int**)malloc(sizeof(int*) * buffer_size); 

    while(true){
        vector<string> precluster = parser.next_cluster(infile);
        if(precluster.size() == 0) break;
        
        vector<int64_t> distribution = processGroup(precluster, 16, kmer_k, true, 2);
        vector<int64_t> rank = computeRankVector(distribution);
        
        if(nPreclusters == buffer_size){
            // Double buffer
            buffer_size *= 2;
            distributions = (int**)realloc(distributions, buffer_size*sizeof(int*));
            ranks = (int**)realloc(ranks, buffer_size*sizeof(int*));
            Component = (int**)realloc(Component, buffer_size*sizeof(int*));
        }
        
        distributions[nPreclusters] = (int*)malloc(sizeof(int) * distribution.size());
        ranks[nPreclusters] = (int*)malloc(sizeof(int) * distribution.size());
        Component[nPreclusters] = (int*)calloc(GenoNum,sizeof(int));
        
        for(int64_t j = 0; j < distribution.size(); j++){
            distributions[nPreclusters][j] = distribution[j];
            ranks[nPreclusters][j] = rank[j];
        }
        nPreclusters++;
    }
        
    int* SixTmer = (int*)calloc(nPreclusters, sizeof(int)); // Some memory for metacluster
    
    MetaCluster meta(kmer_k, nPreclusters, ranks, distributions, SixTmer, GenoNum, Component, 0 , 10, 1);
    
    double MC3_Thresh = 0.94;
    meta.iterMeta(10,MC3_Thresh);
    
    // Free Metacluster parameters
    for(int64_t i = 0; i < nPreclusters; i++){
        free(distributions[i]); 
        free(ranks[i]);
        free(Component[i]);
    }
    free(distributions); 
    free(ranks);
    free(Component);
    free(SixTmer);
    
    int nClasses = *max_element(meta.best, meta.best + nPreclusters) + 1;
    
    // Write preclusters
    vector<ofstream*> outputs;
    for(int64_t i = 0; i < nClasses; i++){
        ofstream* out = new ofstream(workspace + "/cluster_" + to_string(i));
        if(!out->good()){
            cerr << "Error opening file " << workspace + "/cluster_" + to_string(i) << endl;
            exit(1);
        }
        outputs.push_back(out);
    }
    
    ifstream precluster_stream(precluster_file);
    for(int64_t i = 0; i < nPreclusters; i++){
        int64_t cluster = meta.best[i];
        for(auto& read : parser.next_cluster(precluster_stream)){
            *(outputs[cluster]) << ">\n" << read << "\n";
        }
    }
    
    // Clean up
    for(int i = 0; i < nClasses; i++){
        outputs[i]->flush();
        delete outputs[i];
    }
        
    // Concatenate all clusters into one file 
    ofstream results(results_file);
    for(int64_t i = 0; i < nClasses; i++){
        ifstream in(workspace + "/cluster_" + to_string(i));
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
        ofstream X(workspace + "/cluster_" + to_string(i));
        X.close();
    }
    

    cerr << getTimeString() << " Number of final clusters: " << nClasses << endl;
    return nClasses;

}