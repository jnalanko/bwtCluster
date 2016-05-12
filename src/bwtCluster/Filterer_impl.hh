#include <thread>
#include <utility>
#include <atomic>
#include <algorithm>
#include "Algorithms.hh"
#include "Filterer.hh"
#include "fasta_tools.hh"

template <class t_bitvector>
Filterer<t_bitvector>::Filterer(BD_BWT_index<t_bitvector>* index, std::vector<BWT_inversion_iterator<t_bitvector> > iterators, string workspace_path, string paired_end_file) :
        index(index), nThreads(iterators.size()), iterators(iterators), workspace_path(workspace_path), paired_end_file(paired_end_file) {
}

int64_t number_of_good_reads_outputted = 0; // TODO: OUT OF GLOBAL SCOPE TODO TODO
int64_t number_of_bad_reads_outputted = 0;

template <class t_bitvector>
void Filterer<t_bitvector>::output_marked_reads(BWT_inversion_iterator<t_bitvector> iter,
                        std::ostream &good_reads, std::ostream &bad_reads, 
                        std::vector<int64_t>& good_ids,
                        const sdsl::bit_vector& marks){
    
    string newRead = "";
    bool keptAtLeastOne = false;
    bool good = false;
    
    // Collect all reads that have been flagged    
    while(iter.next()){
        if(marks[iter.pos] == 1) good = true;
        if(iter.current_char != '$') newRead.push_back(iter.current_char);
        if(iter.current_char == '$'){
            std::reverse(newRead.begin(), newRead.end());
            // Write the current read if its bit is marked
            output_mutex.lock();
            if(good){
                number_of_good_reads_outputted++;
                good_reads << ">\n" << newRead << "\n";
                good_ids.push_back(iter.read_id);
                keptAtLeastOne = true;
            }
            else{
                number_of_bad_reads_outputted++;
                bad_reads << ">\n" << newRead << "\n";
                // bad_ids << iter.read_id << "\n"; //TODO
            }
            output_mutex.unlock();

            newRead = "";
            good = false;
        }        
    }
    
    
    if(!keptAtLeastOne){
        throw runtime_error(string("Error at "  + string(__FILE__) + " line " +
                to_string(__LINE__) + ": Removed all reads from sample"));
    }
}

template <class t_bitvector>
vector<int64_t> Filterer<t_bitvector>::filterReadsRC(int k, int tau,
            std::string good_out, std::string bad_out){
    
    cerr << getTimeString() << " Filtering reads (RC)" << endl;
    cerr << getTimeString() << " Computing the intervals of right-maximal " << k << "-mers" << endl;

    sdsl::bit_vector marks = get_right_maximal_kmers_rc(*index, k, tau, nThreads);
    
    cerr << getTimeString() << " Inverting the BWT" << endl;
    
    ofstream good_stream(good_out);
    ofstream bad_stream(bad_out);
    
    vector<thread> threads(nThreads);
    vector<int64_t> good_ids_vector; // List of ranks of all reads that were kept
    for(int64_t t = 0; t < nThreads; t++){
        thread x(&Filterer::output_marked_reads, this, iterators[t], ref(good_stream), ref(bad_stream), 
                        ref(good_ids_vector), ref(marks));
        threads[t] = move(x);
    }
    
    for(auto& x : threads) x.join();
    
    good_stream.close(); bad_stream.close();
    return good_ids_vector;
}


template <class t_bitvector>
vector<int64_t> Filterer<t_bitvector>::filterReads(int k, int tau,
            std::string good_out, std::string bad_out){
    
    ofstream good_stream(good_out);
    ofstream bad_stream(bad_out);
    
    cerr << getTimeString() << " Filtering reads (" << nThreads << " threads)" << endl;
    cerr << getTimeString() << " Computing the intervals of right-maximal " << k << "-mers" << endl;
    
    sdsl::bit_vector marks = get_right_maximal_kmers(*index, k, tau, true, nThreads);
    cerr << getTimeString() << " Inverting the BWT" << endl;
    
    vector<thread> threads(nThreads);
    vector<int64_t> good_ids_vector; // List of ranks of all reads that were kept
    for(int64_t t = 0; t < nThreads; t++){
        thread x(&Filterer::output_marked_reads, this, iterators[t], ref(good_stream), ref(bad_stream), 
                        ref(good_ids_vector), ref(marks));
        threads[t] = move(x);
    }
    for(auto& x : threads) x.join();
    
    good_stream.close(); bad_stream.close();
    return good_ids_vector;
}