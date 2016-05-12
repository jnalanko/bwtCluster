#include "Parser.hh"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>

using namespace std;

// TODO: more space efficient

int main(int argc, char** argv){
    if(argc != 3){
        cerr << "Usage: ./preprocessor in.fna out.fna" << endl;
        exit(1);
    }
    
    Parser P;
    ifstream input;
    input.open(argv[1]);
    vector<pair<string,string> > reads = P.parseToVectorWithHeaders(input);
    sort(reads.begin(), reads.end());
    
    ofstream output;
    output.open(argv[2]);
    
    vector<pair<string,string> > reads2;
    int64_t nReadsRemoved = 0;
    for(int64_t i = 0; i < reads.size(); i++){
        if(i != 0 && reads[i].first == reads[i-1].first) {
            nReadsRemoved++;
        }
        else{
            reads2.push_back(reads[i]);
        }
    }
    
    output << reads2.back().second << "\n"; // Header
    output << reads2.back().first << "\n"; // Data
    
    for(int64_t i = 0; i < reads2.size() - 1; i++){
        output << reads2[i].second << "\n"; // Header
        output << reads2[i].first << "\n"; // Data
    }
    
    cerr << "Removed " << nReadsRemoved << " reads" << endl;
}
