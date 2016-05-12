#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include "Parser.hh"

using namespace std;

int main(int argc, char** argv){
    if(argc != 2) return 1;
    
    Parser p;
    ifstream input(argv[1]);
    while(true){
        pair<string,string> x = p.next_read(input);
        if(x.first == "" && x.second == "") break;
        cout << x.second << "\n" << x.first << "\n";
        cout << "---\n";
    }
}