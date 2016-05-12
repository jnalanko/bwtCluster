#include "Parser.hh"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv){
    if(argc != 2 && argc != 3){
        cerr << "Trims 'N' and 'n' characeters off the start and end of all reads" << endl;
        cerr << "Usage: ./trimmer in.fna out.fna" << endl;
        return 1;
    }
    
    Parser parser;
    
    ifstream input(argv[1]);
    if(!input.good()){
        cerr << "Error opening " << argv[1] << endl;
        return 1;
    }
    
    ofstream output(argv[2]);
    if(!output.good()){
        cerr << "Error opening " << argv[2] << endl;
        return 1;
    }
    
    while(true){
        pair<string,string> x = parser.next_read(input);
        if(x.first == "" && x.second == "") break;
        string read = x.first;
        string header = x.second;
        int left = 0; // Leftmost non-N index
        int right = read.size()-1; // Rightmost non-N index
        while(left < read.size() && (read[left] == 'N' || read[left] == 'n'))
            left++;
        while(right >= 0 && (read[right] == 'N' || read[right] == 'n'))
            right--;
        if(left <= right){
            output << header << "\n";
            output << read.substr(left, right-left+1) << "\n";
        } else{
            cerr << "Read discarded:\n" << header << "\n" << read << "\n";
        }
    }
    
}