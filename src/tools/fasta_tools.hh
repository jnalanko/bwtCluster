#ifndef FASTA_TOOLS
#define FASTA_TOOLS

#include <fstream>
#include <utility>
#include <string>
#include <algorithm>
#include "Parser.hh"
#include "tools.hh"

/** @brief Sorts the given fasta file
 * 
 * Sorts the reads in the given fasta file lexicographically. Replaces the content
 * of the headers of the reads with the rank of the read in the original file (zero based).
 * 
 * @param file_in  The filename of the input fasta file
 * @param file_out The filename of the sorted fasta file. Can be the same as file_in
 * 
 * */
void sort_fasta(std::string file_in, std::string file_out);

/** @brief Sorts the given fasta file, and removes duplicate reads
 * 
 * Sorts the reads in the given fasta file lexicographically, and
 * if there are reads with an indentical sequence of nucliotides, keeps only one read
 * from each such set. Replaces the content of the headers of the reads with the rank of 
 * the read in the original file (zero based)
 * 
 * @param file_in  The filename of the input fasta file
 * @param file_out The filename of the output fasta file. Can be the same as file_in
 * 
 * */
void sort_remove_duplicates_update_mate_pairs(std::string fasta_in, std::string fasta_out, std::string paired_in, std::string paired_out);

/** @brief Modifies the read set so that ropebwt outputs the correct forward-BWT
 * 
 * Sorts the read set and rotates the set cyclically forward by one read. Optionally
 * reverses the nucliotide sequences of all reads.
 * 
 * @param input_fasta_filename  The filename of the input fasta file
 * @param output_fasta_filename The filename of the output fasta file. Can be the same as input_fasta_filename
 * @param reverse               True if and only if the read set should be reversed
 * 
 * */
void ropebwt_preprocess(std::string input_fasta_filename, std::string output_fasta_filename, bool reverse);

/** @brief Reverse-complements the fasta file
 * 
 * Computes the reverse-complement of each read in the data.
 * 
 * @param in  The input stream in fasta format
 * @param out The output stream in fasta format
 * 
 * */
void rc_fasta(std::istream& in, std::ostream& out);

/** @brief Rotates the read set forward cyclically
 * 
 * Moves the last *amount* read to the start of the fasta file
 * 
 * @param in  The input stream in fasta format
 * @param out The output stream in fasta format
 * @param amount   The number of reads to shift by
 * 
 * */
void rotate_fasta(std::istream& in, std::ostream& out, int64_t amount);

/** @brief Reverses all read in the input
 * 
 * Outputs the reverse string of each read (note: does not complement the bases)
 * 
 * @param in  The input stream in fasta format
 * @param out The output stream in fasta format
 * 
 * */
void reverse_fasta(std::string input_fasta_filename, std::string output_fasta_filename);

void collect_clusters(std::string input_fasta_filename, std::string clusters_filename, std::string output_fasta_filename);

#endif