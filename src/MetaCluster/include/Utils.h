#ifndef __UTILS_H_

#define __UTILS_H_

#include <cstdio>
#include <vector>
#include <string>

enum OptionType { OptionBool, OptionInt, OptionDouble, OptionString };
enum ParameterType { SIMPLE, INTEGER, FLOAT, STRING, };

void AddOption(const std::string &long_name, const std::string &short_name,
        bool &bool_option, const std::string &description);
void AddOption(const std::string &long_name, const std::string &short_name,
        int &int_option, const std::string &description);
void AddOption(const std::string &long_name, const std::string &short_name,
        double &double_option, const std::string &description);
void AddOption(const std::string &long_name, const std::string &short_name,
        std::string &string_option, const std::string &description);

void ProcessOptions(int &argc, char *argv[]);
std::string OptionDescriptions();

void AddParameter(const char *name, void *pointer, ParameterType type);
void ProcessParameters(int &argc, char *argv[]);

#endif
