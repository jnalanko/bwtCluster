#include "tools.hh"
#include <string>
#include <cassert>
#include <ctime>
#include <chrono>
#include <fstream>
#include <cstdio>
#include <cerrno>

using namespace std;

int64_t binaryToInt(string binary){
    // Assume most significant bit first

    int64_t result = 0;
    assert(binary.size() <= 64);
    for(int i = 0; i < binary.size(); i++){
        char c = binary[binary.size() - 1 - i];
        assert(c == '0' || c == '1');
        if(c == '1') result += ((int64_t)1 << i);
    }
    return result;
}


string getTimeString(){
    std::time_t result = std::time(NULL);
    string time = std::asctime(std::localtime(&result));
    return time.substr(0,time.size() - 1); // Trim the trailing newline
}

int64_t getUnixEpoch(){
    auto time = std::chrono::system_clock::now().time_since_epoch();
    return std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
}

void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

void trim(std::string &s) {
    ltrim(s); rtrim(s);
}

void write_to_disk(string& text, string filepath){
    ofstream stream(filepath);
    stream << text;
    stream.flush();
}

char* read_from_disk_c_string(std::string filename)
{
  std::FILE *fp = std::fopen(filename.c_str(), "rb");
  if (fp)
  {
    std::fseek(fp, 0, SEEK_END);
    int64_t size = std::ftell(fp);
    char* contents = (char*)malloc(size + 1);
    std::rewind(fp);
    int unused = std::fread(&contents[0], 1, size, fp);
    (void) unused; // TODO: do something with this value
    std::fclose(fp);
    return contents;
  }
  throw(errno);
}

std::string read_from_disk(std::string filename)
{
  std::FILE *fp = std::fopen(filename.c_str(), "rb");
  if (fp)
  {
    std::string contents;
    std::fseek(fp, 0, SEEK_END);
    contents.resize(std::ftell(fp));
    std::rewind(fp);
    int unused = std::fread(&contents[0], 1, contents.size(), fp);
    (void) unused; // TODO: do something with this value
    std::fclose(fp);
    return contents;
  }
  throw(errno);
}

void copy_file(std::string file_from, std::string file_to){
    ifstream from(file_from);
    ofstream to(file_to);
    string line;
    while(getline(from, line)) 
        to << line << "\n";
}