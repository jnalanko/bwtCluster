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
#include "TomAlgorithm.h"

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

int main(int argc, char** argv){
    if(argc != 3){
        cerr << "Usage: ./program [precluster_file] [k]" << endl;
        return 1;
    }
    string precluster_file = argv[1];
    int64_t kmer_k = stoi(argv[2]);
    
    Parser parser;
    ifstream infile(precluster_file);
    vector<vector<string> > preclusters = parser.parse_preclusters(infile);

    int GenoNum = 1;
    int** distributions = (int**)malloc(sizeof(int*) * preclusters.size()); // todo: int64_t. Is now int because Metacluster takes int**
    int** ranks = (int**)malloc(sizeof(int*) * preclusters.size()); // todo: int64_t
    int** Component = (int**)malloc(sizeof(int*) * preclusters.size());
    
    cerr << "Distributions take: " << preclusters.size() * ((kmer_k*2) << 1) * sizeof(int) << endl;
    int* ReverPosition = getReverPosition(kmer_k);
    int ReverSize = tomReverSize(kmer_k);
    int NoReverSize = 1<<(kmer_k<<1);
            
    //#pragma omp parallel for 
    for(int64_t i = 0; i < preclusters.size(); i++){
        cout << i << endl;
        vector<int64_t> distribution = processGroup(preclusters[i], 16, kmer_k, false, 2);
        ranks[i] = toSpear(tomNormalize_rever(&(distribution[0]),ReverPosition,kmer_k,false),ReverSize);
        distributions[i] = (int*)malloc(sizeof(int) * distribution.size());
         for(int j=0;j<NoReverSize;++j)
            distributions[i][j] = distribution[j];
        for(int j=0;j<NoReverSize;++j)
        {
            int rever = tomReverComple(j,kmer_k);
            if(j<=rever)
            {
                distributions[i][j] += distributions[i][rever];
                distributions[i][rever] = distributions[i][j];
            }
        }
        Component[i] = (int*)calloc(GenoNum,sizeof(int));
        /*
        vector<int64_t> rank = computeRankVector(distribution);
        
        distributions[i] = (int*)malloc(sizeof(int) * distribution.size());
        ranks[i] = (int*)malloc(sizeof(int) * distribution.size());
        Component[i] = (int*)calloc(GenoNum,sizeof(int));
        for(int64_t j = 0; j < distribution.size(); j++){
            distributions[i][j] = distribution[j];
            ranks[i][j] = rank[j];
        }*/
    }
    
    int* SixTmer = (int*)calloc(preclusters.size(), sizeof(int));
    
    MetaCluster meta(kmer_k, preclusters.size(), ranks, distributions, SixTmer, GenoNum, Component, 0 , 6, 1);
    
    double MC3_Thresh = 0.94;
    meta.iterMeta(10,MC3_Thresh);
    
    /*for(int i = 0; i < preclusters.size(); i++){
        cout << meta.best[i] << endl;
    }*/
    int nClasses = *max_element(meta.best, meta.best + preclusters.size()) + 1;
    cout << "Classes in the end: " << nClasses << endl;
    
    vector<ofstream*> outputs;
    for(int64_t i = 0; i < nClasses; i++){
        ofstream* out = new ofstream("./mc_out/cluster_" + to_string(i));
        if(!out->good()){
            cerr << "Error opening file " << "./mc_out/cluster_" + to_string(i) << endl;
            exit(1);
        }
        outputs.push_back(out);
    }
    
    for(int64_t i = 0; i < preclusters.size(); i++){
        int64_t cluster = meta.best[i];
        for(auto& read : preclusters[i]){
            *(outputs[cluster]) << ">\n" << read << "\n";
        }
    }
    
    // Clean up
    for(int i = 0; i < nClasses; i++){
        outputs[i]->flush();
        delete outputs[i];
    }
    
    /*
    for(int64_t i = 0; i < preclusters.size(); i++){
        free(distributions[i]); 
        free(ranks[i]);
        free(Component[i]);
    }
    free(distributions); 
    free(ranks);
    free(Component);
    free(SixTmer); */ //TODO
    //delete[] ReverPosition; // Do we have to do this?
    
//explicit MetaCluster(int KmerLen_,int size_,int**spear_,int**KmerDistri_,int*SixTmer_,int GenoNum_,int** Component_,int kmeansize,int MaxSpecies_,int MinSpecies_)
}

/*
int KmerLen_; // K-mer length
int size_; // Number of preclusters?
int** spear_; // Array of rank vectors (size x 4^KmerLen)?
int** KmerDistri_; // Array of Kmer distribution vectors (size x 4^KmerLen)?
int* SixTmer_; // Array of length size, irrelevant in the end?
int GenoNum_; // Always 1??
int** Component_; (size x GenoNum) array of all zeroes
int kmeansize; // Initial k-means k. Will be iterated from here. If 0, use MaxSpecies.
int MaxSpecies_; // Upper bound for kmeans k
int MinSpecies; // Lower bound for kmeans k
*/
