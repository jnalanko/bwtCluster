#include "bwt.hh"
#include <cstring>
#include <iostream>
#include "tools.hh"
#include <fstream>
#include <string>

extern "C" { 
#include "dbwt.h"
}

extern "C" { 
#include "ropebwt2.h"
}

using namespace std;

char* bwt_dbwt(char* text){

    unsigned int last;
    long n = strlen(text);
    char* d = (char*)dbwt_bwt((uchar*)text, n, &last, 0);
    d = (char*)realloc(d, (n + 2) * sizeof(char));
    d[last] = '$';
    d[n + 1] = 0;

    return d;
}

void bwt_ropebwt(std::string input_filename, string bwt_path, std::string work_directory, int64_t memory_bytes){
    cerr << getTimeString() << " Running ropebwt" << endl;
    string temp = work_directory + "/temp.txt";
    int argc = 7;
    char* argv[7];
    for(int i = 0; i < argc; i++) argv[i] = (char*)malloc(sizeof(char)*100);
    
    strcpy(argv[0], "ropebwt2");
    strcpy(argv[1], input_filename.c_str());
    strcpy(argv[2], "-R");
    strcpy(argv[3], "-m");
    strcpy(argv[4], to_string(memory_bytes).c_str());
    strcpy(argv[5], "-o");
    strcpy(argv[6], temp.c_str());
    main_ropebwt2(argc,argv);
    for(int i = 0; i < argc; i++) free(argv[i]);
    /*if(system((ropebwt_path + " " + input_filename + " -R -m " + to_string(extra_memory_megabytes) + " -o " + temp + " 2> /dev/null").c_str())){
        cerr << "Error executing: " << ropebwt_path + " " + input_filename + " -R -m " + to_string(extra_memory_megabytes) + " -o " + temp + " 2> /dev/null" << endl;
        exit(1);
    }*/
    
    // To match old psascan based bwt, add a dollar to the second position
    // Also remove the endline in the end
    ifstream in(temp);
    string s; in >> s; // todo: memory overhead?
    ofstream out(bwt_path);
    out << s[0] << '$'; 
    for(int64_t i = 1; i < s.size(); i++)
        out << s[i];
    out.flush();
    //system(("perl -pi -e 'chomp if eof' " + bwt_path).c_str()); // Remove endline from the end of the file
}


string bwt_pSAscan(string& input, string work_directory, string pSAscan_path, string bwt_path, int64_t memory_budget_megabytes){

    // Clean up
    if(system(("rm " + work_directory + "/string.txt 2> /dev/null").c_str())){
        cerr << "Error executing: " << "rm " + work_directory + "/string.txt 2> /dev/null" << endl;
        exit(1);
    }
    if(system(("rm " + work_directory + "/string.txt.sa5 2> /dev/null").c_str())){
        cerr << "Error executing: " << "rm " + work_directory + "/string.txt.sa5 2> /dev/null" << endl;
        exit(1);        
    }
    
    int64_t n = input.size();

    // Write the input into a file for SAscan
    input.push_back('$');
    write_to_disk(input, work_directory + "/string.txt");

    // Build the suffix array giving n bytes of memory to sascan
    cerr << getTimeString() << " Running parallel sascan" << endl;
    if(!system((pSAscan_path + " " + work_directory + "/string.txt --mem=" + to_string(memory_budget_megabytes) + " 2> /dev/null").c_str())){
        cerr << "Error executing: " << "rm " + work_directory + "/string.txt.sa5 2> /dev/null" << endl;
        exit(1);
    }

    // Build the BWT from the SA. The SA elements are represented
    // by 5 bytes each in the file
    FILE* fptr = fopen((work_directory + "/string.txt.sa5").c_str(),"rb");
    if(fptr==NULL){
      printf("IO Error at bwt_SA_scan");
      exit(1);
    }

    int64_t x = 0;

    cerr << getTimeString() << " Streaming suffix array from disk" << endl;
    ofstream bwt_stream(bwt_path);
    for(int64_t i = 0; i < (n+1) * 5; i++){
        unsigned char c = (unsigned char)fgetc(fptr);
        x |= ((int64_t)c) << (8 * (i%5)); // Shift the char into position
        if((i + 1) % 5 == 0){
            int64_t SA_value = x;
            bwt_stream << (input[SA_value == 0 ? n : SA_value - 1]);
            x = 0;
        }
    }
    bwt_stream.flush();
    input.pop_back(); // The dollar
    return bwt_path;
}
