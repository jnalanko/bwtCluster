#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>
#include "Parser.hh"
#include "tools.hh"

using namespace std;

// http://www.cplusplus.com/forum/general/10898/
std::string dec2bin(int64_t n){
    const int size=sizeof(n)*8;
    std::string res;
    bool s = 0;
    for (int a = 0; a < size; a++){
        bool bit = n >> (size-1);
        if (bit)
            s = 1;
        if (s)
            res.push_back(bit + '0');
        n <<= 1;
    }
    if (!res.size())
        res.push_back('0');
    return res;
}

struct Parser::parseConcatenateWithHeaders_return Parser::parseConcatenateWithHeaders(std::istream& input){
    struct parseConcatenateWithHeaders_return result;
    while(true){
        string line; getline(input,line); trim(line);
        if(line.size() == 0 && !input.eof()) continue;
        if(input.eof()) break;
        if(line[0] == '>'){
            if(result.starts.size() != 0){
                int64_t read_end = result.data.size() - 1;
                int64_t read_start = result.starts.back() + result.headerLengths.back();
                cout << "start: " << read_start << ", end: " << read_end << endl;
                result.readLengths.push_back(read_end - read_start + 1);
            }
            result.starts.push_back(result.data.size());
            result.headerLengths.push_back(line.size());
            
        }
        result.data += line;
    }
    int64_t last_read_end = result.data.size() - 1;
    int64_t last_read_start = result.starts.back() + result.headerLengths.back();
    result.readLengths.push_back(last_read_end - last_read_start + 1);
    result.data.shrink_to_fit();
    return result;
}

pair<string,int64_t> Parser::parseConcatenate(std::istream &input, char separator, string paired_end_file){
    ofstream read_pairs;
    if(paired_end_file != "") read_pairs.open(paired_end_file);
    
    string data = "";
    int64_t nReads = 0;
    string read = "";
    vector<char> alphabet = {'A','C','G','T'};
    bool readGood = true; // Does the read contain only letters from the allowed alphabet?
    
    int64_t number_parsed = 0; // Number of reads parsed, regardless whether hey had a valid alphabet
    bool previous_good = false; // Did the previous read contain only letters from the allowed alphabet
    while(true){
        
        string line; getline(input,line); trim(line);
        if(line.size() == 0 && !input.eof()) continue;
        if(input.eof() || line[0] == '>'){ // Add the current read to the data if it's good
            if(read.size() != 0) number_parsed++; // read.size() == 0: Happens at least when reading the first line of the file
            if(read.size() != 0 && readGood){                
                data += separator + read;
                nReads++;
                
                // Write the mate pair if appropriate
                if(number_parsed % 2 == 0 && previous_good && paired_end_file != ""){
                    read_pairs << nReads - 2 << " " << nReads - 1 << "\n";
                }
                
                previous_good = true;
            } else previous_good = false;
         
            read = "";
            readGood = true;
            
            if(input.eof()){
                if(number_parsed % 2 != 0 && paired_end_file != ""){
                    cerr << "Error: odd number of reads in input, but paired ends enabled" << endl;
                    exit(1);
                }
                data.shrink_to_fit();
                return make_pair(data,nReads);
            }
            
        } else{ // Add the characters on the line to the current read
            for(char c : line){
                c = toupper(c);
                if(find(alphabet.begin(), alphabet.end(), c) == alphabet.end())
                    readGood = false;
                read.push_back(c);
            }
        }
    }
}

std::vector<std::string> Parser::parseToVector(std::istream &input){
    vector<string> v;
    string read = "";
    vector<char> alphabet = {'A','C','G','T'};
    bool readGood = true; // Does the read contain only letters from the allowed alphabet?
    while(true){
        string line; getline(input,line); trim(line);
        if(line.size() == 0 && !input.eof()) continue;
        if(input.eof() || line[0] == '>'){ // Add the current read to the vector if it's good
            if(read.size() != 0 && readGood) {
                read.shrink_to_fit();
                v.push_back(read);
            }
                
            read = "";
            readGood = true;
            
            if(input.eof()) return v;
            
        } else{ // Add the characters on the line to the current read
            for(char c : line){
                c = toupper(c);
                if(find(alphabet.begin(), alphabet.end(), c) == alphabet.end())
                    readGood = false;
                read.push_back(c);
            }
        }
    }
}

std::vector<std::pair<std::string, std::string> > Parser::parseToVectorWithHeaders(std::istream &input){
    
    vector<char> alphabet = {'A','C','G','T'};
    vector<string> lines;
    
    // Read the non-empty lines
    while(true){
        string line; getline(input,line); trim(line);
        if(input.eof()) break;
        if(line.size() != 0)
        lines.push_back(line);
    }
    
    vector<pair<string,string> > result;
    string header = "";
    string read = "";
    for(int64_t i = 0; i < lines.size(); i++){
        string line = lines[i];
        if(line[0] == '>') header = line;
        else{
            read += line;
            if(i == lines.size() - 1 || lines[i+1][0] == '>'){
                bool read_has_valid_chars = true;
                for(int64_t j = 0; j < read.size(); j++){
                    read[j] = toupper(read[j]);
                    if(find(alphabet.begin(), alphabet.end(), read[j]) == alphabet.end())
                        read_has_valid_chars = false;
                }
                if(read_has_valid_chars){
                    read.shrink_to_fit();
                    header.shrink_to_fit();
                    result.push_back({read,header});
                }
                header = ""; read = "";
            }
        }
    }
    
    return result;  
}


std::vector<std::vector<int64_t> > Parser::parseGroups(std::istream &groups){
    vector<vector<int64_t> > result;
    while(true){
        if(groups.eof()) break;
        string line;
        getline(groups,line);
        if(line.size() == 0) continue;
        stringstream lineStream(line);
        int64_t x;
        vector<int64_t> v;
        while(lineStream >> x){
            v.push_back(x);
        }
        result.push_back(v);
    }
    return result;
}

std::vector<std::vector<string> > Parser::parse_clusters(std::istream &input){
    vector<pair<string, string> > asd = parseToVectorWithHeaders(input); // (read, header) pairs
    vector<vector<string> > preclusters;
    for(int64_t i = 0; i < asd.size(); i++){
        string read = asd[i].first;
        string header = asd[i].second;
        if(header == ">cluster_start"){
            vector<string> v;
            preclusters.push_back(v);
        }
        preclusters.back().push_back(read);
    }
    return preclusters;
}

std::vector<std::string> Parser::next_cluster(std::istream &input){
    vector<string> cluster;
    if(input.eof()){
        return cluster;
    }
    string read = "";
    while(true){
        string line; getline(input,line); trim(line);
        if(input.eof() || (line.substr(0,14) == ">cluster_start" && cluster.size() != 0)){
            cluster.push_back(read);
            return cluster;
        }
        
        if(line[0] == '>'){
            if(read != "") cluster.push_back(read); // read == "" if this is the first read of the file
            read = "";
        } else{
            read += line;
        }
    }
}

std::pair<std::string, std::string> Parser::next_read(std::istream &input){
    string header = "";
    string read = "";
    getline(input,header);
    if(!input.good()) return {"", ""};
    while(input.peek() != '>'){
        string buffer;
        getline(input,buffer);
        if(!input.good()) break;
        read += buffer;
    }
    return {read, header};
}


