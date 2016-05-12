#ifndef PARSER_HH
#define PARSER_HH

#include <iostream>
#include <string>
#include <vector>
#include <utility>

class Parser{
  
private:
    
public:


    /* @brief    Function parseConcatenate
     *
     * Parses the reads in the input stream by adding a separator character to the start of
     * each read and concatenating the reads. Makes all characters upper case. 
     * Discards reads that have characters outside of the alphabet {A,C,G,T}.
     * If paired_end_file is not the empty string, Assumes that every two consecutive reads in the file are mate pairs.
     * If this is the case, Writes the ranks of the mate pairs of the reads in the concatenation (i.e. the reads that had a valid alphabet)
     * into the file paired_end_file such that on each line of the file there is a pair of integers (seperated by
     * a space), that specifies the ranks of the reads of a mate pair in the concatenation.
     *
     * @param    input           An input stream to FASTA format read data
     * @param    separator       The character to add to the start of each read
     * @param    paired_end_file A file to write the mate pairs
     * @return                   A pair with the concatenated data and the number of reads in the data
     */
    std::pair<std::string, int64_t> parseConcatenate(std::istream &input, char separator, std::string paired_end_file = "");    
    
    /* @brief    Function parseToVector
     *
     * Parses the reads in the input into a std::vector. Makes all characters upper case.
     * Discards reads that have characters outside of the alphabet {A,C,G,T}.
     *
     * @param    input        An input stream to FASTA format read data
     * @return                A vector with all the reads
     */
    std::vector<std::string> parseToVector(std::istream &input);

    /* @brief    Function parseToVector
     *
     * Parses the reads in the input into a std::vector of pairs of strings, the first being
     * the DNA sequence, and the second the fasta header. Makes all characters upper case in the DNA.
     * Discards reads that have characters outside of the alphabet {A,C,G,T}.
     *
     * @param    input        An input stream to FASTA format read data
     * @return                A vector with all the reads
     */
    std::vector<std::pair<std::string, std::string> > parseToVectorWithHeaders(std::istream &input);
    
    std::vector<std::vector<std::string> > parse_clusters(std::istream &input);

    /* @brief    Function parseGroups
     *
     * Parses an input stream format that contains one set of space-separated integers for each line
     *
     * @param    input        An input stream to the data
     * @return                A vector containing the set integers for each line
     */
    std::vector<std::vector<int64_t> > parseGroups(std::istream &groups);

    std::vector<std::string> next_cluster(std::istream &input);
    
    std::pair<std::string, std::string> next_read(std::istream &input);
    
    int64_t count_clusters(std::istream& input);
    
    struct parseConcatenateWithHeaders_return{
        std::string data; // header, data, header, data...
        std::vector<int64_t> starts;
        std::vector<int64_t> headerLengths;
        std::vector<int64_t> readLengths;
    };
    
    struct parseConcatenateWithHeaders_return parseConcatenateWithHeaders(std::istream& input);


};

#endif
