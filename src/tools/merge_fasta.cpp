#include <iostream>
#include <string>
#include "Parser.hh"
#include <utility>
#include <fstream>
#include <stdint.h>


using namespace std;

vector<pair<string,string> > read_file(string filename){

    ifstream input(filename);
    if(!input.good()){
        cerr << "Error opening file " << filename << endl;
        exit(-1);
    }

    string header, read;
    vector<pair<string, string> > v;
    int state = 0;

    string line;
    while(true){
        getline(input,line);
        if(!input.good()) break;
        if(state == 0){
            header = line;
            state = 1;
        } else{
            read = line;
            state = 0;
            v.push_back({header,line});
        }
    }
    return v;
}

int main(int argc, char** argv){
    if(argc != 4){
        cerr << "Usage: merge_fasta file1.fna file2.fna out.fna" << endl;
        cerr << "Merges file1.fna and file2.fna into out.fna taking reads alternatively" << endl;
        cerr << "from both files. ASSUMES ALL READS TAKE ONLY A SINGLE LINE EACH!!!" << endl;
        return -1;
    }
    
    vector<pair<string, string> > reads1, reads2; // (header, read) pairs
    reads1 = read_file(argv[1]);
    reads2 = read_file(argv[2]);
    ofstream output(argv[3]);
    if(!output.good()){
        cerr << "Error opening file " << argv[3] << endl;
        return -1;
    }

    int64_t p1 = 0; // index in reads1
    int64_t p2 = 0; // index in reads2

    int state = 0;
    while(p1 < reads1.size() || p2 < reads2.size()){
        if(state == 0){
            output << reads1[p1].first + "/1" << "\n" << reads1[p1].second << "\n";
            p1++;
            if(p2 < reads2.size()) state = 1;
        }
        else if(state == 1){
            output << reads2[p2].first + "/2" << "\n" << reads2[p2].second << "\n";
            p2++;
            if(p1 < reads1.size()) state = 0;
        }
    }
    
    return 0;
}



