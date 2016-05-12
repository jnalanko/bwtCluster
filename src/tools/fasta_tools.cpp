#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <iterator>
#include "tools.hh"
#include "Parser.hh"

using namespace std;

struct Read{
    int64_t start;
    int64_t original_rank;
    int length;
};

// A functor to compare two reads lexicographically
struct compare_reads{
    
    string* data;
    
    compare_reads(string& data) : data(&data) {}
    
    bool operator()(const Read& r1, const Read& r2){
        for(int64_t i = 0; i < min(r1.length, r2.length); i++){
            if((*data)[r1.start+i] < (*data)[r2.start+i]) return true;
            if((*data)[r1.start+i] > (*data)[r2.start+i]) return false;
        }
        
        // The shorter string is a prefix of the longest string
        if(r1.length < r2.length) return true;
        else return false;
    }
};


void do_merge(vector<Read>&reads, vector<Read>& aux, 
              int64_t left_start, int64_t left_end, int64_t right_start, int64_t right_end, 
              compare_reads& smaller_than){
    int64_t left_ptr = left_start;
    int64_t right_ptr = right_start;
    int64_t aux_ptr = left_ptr;
    while(left_ptr < left_end || right_ptr < right_end){
        if(right_ptr >= right_end){
            aux[aux_ptr] = reads[left_ptr];
            left_ptr++;
        }
        else if(left_ptr >= left_end){
            aux[aux_ptr] = reads[right_ptr];
            right_ptr++;
        }
        else{
            if(smaller_than(reads[left_ptr], reads[right_ptr])){
                aux[aux_ptr] = reads[left_ptr];
                left_ptr++;
            }
            else{
                aux[aux_ptr] = reads[right_ptr];
                right_ptr++;
            }
        }
        aux_ptr++;
    }
}

void merge_sort(vector<Read>& reads, compare_reads& smaller_than){
    int64_t n = reads.size();
    vector<Read> aux(n);
    for(int64_t block_size = 1; block_size < n; block_size *= 2){
        
        for(int64_t left_start = 0; left_start < n; left_start += block_size*2){
            int64_t left_end = min(left_start + block_size, n); // Exclusive
            int64_t right_start = left_end;
            int64_t right_end = min(right_start + block_size, n); // Exclusive
            
            do_merge(reads, aux, left_start, left_end, right_start, right_end, smaller_than);
        }

        for(int64_t i = 0; i < aux.size(); i++)
            reads[i] = aux[i];
    }
}

void write_unique_reads(vector<Read>& reads, string& data, string filename){
    ofstream out(filename);
    for(int64_t i = 0; i < reads.size(); i++){
        string r = data.substr(reads[i].start, reads[i].length);
        if(i == 0 || r != data.substr(reads[i-1].start, reads[i-1].length)){
            out << "> rank in original file: " << reads[i].original_rank << "\n" << r << "\n";    
        }
    }
}

void remove_duplicates_in_sorted_list(vector<Read>& reads, string& data){
    int64_t pos = 0;
    for(int64_t i = 0; i < reads.size(); i++){
        string r = data.substr(reads[i].start, reads[i].length);
        if(i == 0 || r != data.substr(reads[i-1].start, reads[i-1].length)){
            reads[pos] = reads[i];
            pos++;
            //out << "> rank in original file: " << reads[i].original_rank << "\n" << r << "\n";    
        }
    }
    reads.resize(pos);
}

vector<Read> build_read_vector(string data){
    vector<Read> reads;
    Read current;
    int64_t rank = 0;
    for(int64_t i = 0; i < data.size(); i++){
        if(data[i] == '$') current.start = i + 1;
        if(data[i+1] == '$' || i == data.size() - 1){
            current.length = i - current.start + 1;
            current.original_rank = rank;
            reads.push_back(current);
            rank++;
        }
    }
    reads.shrink_to_fit();
    return reads;
}

void sort_fasta(std::string file_in, std::string file_out){
    Parser P;
    ifstream in(file_in);
    pair<string, int64_t> parse_result = P.parseConcatenate(in, '$');
    in.close();
    
    string data = parse_result.first;
    vector<Read> reads = build_read_vector(data);
    
    compare_reads smaller_than(data);
    merge_sort(reads, smaller_than);
    
    ofstream out(file_out);
    for(auto r : reads){
        out << ">\n" << data.substr(r.start,r.length) << "\n";
    }
    
    data = ""; data.shrink_to_fit();
}

void update_mate_pairs(vector<int64_t>& mapping, string paired_in_file, string paired_out_file){ 
    
    assert(paired_in_file != paired_out_file);
    
    ifstream pairs_in(paired_in_file);
    ofstream pairs_out(paired_out_file);
    
    while(true){
        string line;
        getline(pairs_in, line);
        if(!pairs_in.good()) break;
        stringstream ss(line);
        int64_t r1, r2;
        ss >> r1 >> r2;
        if(mapping[r1] != -1 && mapping[r2] != -1)
            pairs_out << mapping[r1] << " " << mapping[r2] << "\n";
    }
}


void sort_remove_duplicates_update_mate_pairs(std::string fasta_in, std::string fasta_out, std::string paired_in, std::string paired_out){
    Parser P;
    ifstream in(fasta_in);
    pair<string, int64_t> parse_result = P.parseConcatenate(in, '$');
    in.close();
    
    string data = parse_result.first;
    vector<Read> reads = build_read_vector(data);
    
    compare_reads smaller_than(data);
    merge_sort(reads, smaller_than);
    
    vector<int64_t> perm(reads.size(),-1); // The inverse of the permutation the sorting does on the reads
    remove_duplicates_in_sorted_list(reads, data);
    
    for(int64_t i = 0; i < reads.size(); i++){
        perm[reads[i].original_rank] = i;
    }
    
    if(paired_in != "") update_mate_pairs(perm, paired_in, paired_out);
    
    ofstream out(fasta_out);
    for(auto r : reads){
        out << ">\n" << data.substr(r.start,r.length) << "\n";
    }    
    
    data = ""; data.shrink_to_fit();    
}


void ropebwt_preprocess(std::string input_fasta_filename, std::string output_fasta_filename, bool reverse){
    Parser P;
    ifstream in(input_fasta_filename);
    string data = P.parseConcatenate(in,'$').first;
    data.push_back('$');    
    
    if(reverse) std::reverse(data.begin(), data.end());
    data.pop_back();
    data.shrink_to_fit();
    
    // Collect the read start positions and lengths in data
    vector<Read> reads;
    Read current;
    for(int64_t i = 0; i < data.size(); i++){
        if(data[i] == '$') current.start = i + 1;
        if(data[i+1] == '$' || i == data.size() - 1){
            current.length = i - current.start + 1;
            reads.push_back(current);
        }
    }
    
    compare_reads smaller_than(data);
    merge_sort(reads, smaller_than);
    
    ofstream out(output_fasta_filename);
    
    // Write the last read first
    out << ">\n" << data.substr(reads.back().start, reads.back().length) << "\n";
    
    // Write the rest of the reads
    for(int64_t i = 0; i < reads.size() - 1; i++){
        out << ">\n" << data.substr(reads[i].start, reads[i].length) << "\n";
    }
    out.flush();  
}

// Reverse-complements all reads in the input file
void rc_fasta(std::istream& in, std::ostream& out){
    Parser P;
    vector<pair<string,string> > reads = P.parseToVectorWithHeaders(in);

    for(auto& read : reads){
        out << read.second << " RC \n"; // Header
        
        for(int64_t i = read.first.size() - 1; i >= 0; i--){
            char c = read.first[i];
            if(c == 'A') out << 'T';
            else if(c == 'C') out << 'G';
            else if(c == 'G') out << 'C';
            else if(c == 'T') out << 'A';
            else{
                cerr << "Invalid character: " << c << endl;
                exit(1);
            }
        }
        out << "\n";
    }
    out.flush();
}

// Moves the last *amount* reads to the start of the file
void rotate_fasta(std::istream& in, std::ostream& out, int64_t amount){
    Parser P;
    vector<pair<string,string> > reads = P.parseToVectorWithHeaders(in);
    for(int64_t i = reads.size() - amount; i < reads.size(); i++){
        out << reads[i].second << "\n"; // Header
        out << reads[i].first << "\n"; // Data
    }
    for(int64_t i = 0; i < reads.size() - amount; i++){
        out << reads[i].second << "\n"; // Header
        out << reads[i].first << "\n"; // Data
    }
    out.flush();
}

void reverse_fasta(std::string input_fasta_filename, std::string output_fasta_filename){
    Parser P;
    ifstream in(input_fasta_filename);
    vector<string> data = P.parseToVector(in);
    ofstream out(output_fasta_filename);
    for(auto read : data){
        reverse(read.begin(), read.end());
        out << ">\n" << read << "\n";
    }
    out.flush();
}

vector<string> split(string s){
    istringstream iss(s);
    vector<string> tokens{istream_iterator<string>{iss},
                          istream_iterator<string>{}};
    return tokens;
}

void collect_clusters(std::string input_fasta_filename, std::string clusters_filename, std::string output_fasta_filename){

    Parser P;
    ifstream read_in_file(input_fasta_filename);
    if(!read_in_file.good()){
        cerr << "Error reading file: " << input_fasta_filename << endl;
        exit(1);
    }
    vector<string> reads = P.parseToVector(read_in_file);
    
    ifstream ids_in(clusters_filename);
    if(!ids_in.good()){
        cerr << "Error reading file: " << clusters_filename << endl;
        exit(1);
    }
    
    string line;
    ofstream out(output_fasta_filename);
    if(!out.good()){
        cerr << "Error opening file: " << output_fasta_filename << endl;
        exit(1);
    }
    
    while(getline(ids_in,line)){    
        vector<string> tokens = split(line);
        for(int64_t i = 0; i < tokens.size(); i++){
            string token = tokens[i];
            int64_t id = stoi(token);
            if(i == 0) out << ">cluster_start\n";
            else out << ">\n";
            out << reads[id] << "\n";
        }
    }
}
