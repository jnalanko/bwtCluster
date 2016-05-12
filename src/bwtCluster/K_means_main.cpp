#include <K_means.hh>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char** argv){
    if(argc != 5){
        cerr << "Usage: ./program [preclusterfile.fna] [out_directory without trailing slash] [k] [maxrounds]" << endl;
        return 1;
    }
    
    string infile = argv[1];
    string out_dir = argv[2];
    int64_t k = stoi(argv[3]);
    int64_t maxRounds = stoi(argv[4]);
    K_means K(infile,4);
    K.run();
    
}