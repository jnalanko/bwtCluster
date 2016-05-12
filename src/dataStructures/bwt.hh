#ifndef BWT_HH
#define BWT_HH

#include <string>

// Returns the path to a file that contains the bwt of the input
std::string bwt_pSAscan(std::string& input, std::string work_directory, std::string pSAscan_path, std::string bwt_path, int64_t max_memory_megabytes);

// Returns the path to a file that contains the bwt of the input
void bwt_ropebwt(std::string input_filename, std::string bwt_path,  std::string work_directory, int64_t extra_memory_bytes);

char* bwt_dbwt(char* text);

#endif
