#ifndef TOOLS_HH
#define TOOLS_HH

#include <cstdlib>
#include <string>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

int64_t binaryToInt(std::string binary);
std::string getTimeString();
int64_t getUnixEpoch();
std::string getLabel(struct iterator_state* state);


// trim whitespace from start
void ltrim(std::string &s);

// trim whitespace from end
void rtrim(std::string &s);

// trim whitespace from both ends
void trim(std::string &s);

void copy_file(std::string file_from, std::string file_to);

void write_to_disk(std::string& text, std::string filepath);
char* read_from_disk_c_string(std::string filepath);
std::string read_from_disk(std::string filepath);


#endif
