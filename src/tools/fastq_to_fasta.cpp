#include <iostream>
#include <string>

using namespace std;

int main(int argc, char** argv){
    if(argc != 1){
        cerr << "Usage: pipe fastq into the program, prints fasta to stdout" << endl;
        cerr << "Takes every fourth line in the input starting from the second line" << endl;
        cerr << "Gives each read a unique header" << endl;
        return -1;
    }
    
    int state = 0;
    int64_t read_count = 0;
    string line = "";
    while(true){
        getline(cin,line);
        if(!cin.good()) break;
        if(state == 1){
            cout << ">" << read_count << "\n";
            cout << line << "\n";
            read_count++;
        }
        state = (state + 1) % 4;
    }
       
}
