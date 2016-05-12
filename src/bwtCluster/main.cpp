#include "BWT_cluster.hh"
#include "tools.hh"
#include "Parser.hh"
#include "K_means.hh"
#include "fasta_tools.hh"
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <utility>
#include <sstream>
#include <exception>
#include <map>
#include <set>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include "get_rss.hh"
#include <ctime>

using namespace std;

class Config{
public:
    int read_filter_k, read_filter_threshold, grouping_k, max_group_size, k_means_rounds, 
        min_group_size, k_means_k_mer_length, k_means_number_of_clusters, number_of_threads;
    string input_reads, good_reads_out, bad_reads_out, final_clusters_out, preclusters_out,
           workspace_path, bitvector;
    bool paired_ends, run_metacluster, only_kmeans;
    
    bool read_bool(string s){
        if(s != "true" && s != "false"){
            cerr << "Invalid config parameter: " << s << endl;
            exit(1);                   
        }
        return (s == "true" ? true : false);    
    }
    
    // Names of the parameters in the config file and as command line parameters:
    const vector<string> names = {"--read-filtering-k", "--read-filtering-threshold", "--preclustering-k", "--precluster-min-size", "--precluster-max-size", 
        "--kmeans-kmer-length", "--kmeans-number-of-clusters", "--threads",  "--filtering-good-out", "--filtering-bad-out", 
         "--preclusters-out", "--final-clusters-out", "--workspace", "--bitvector", "--paired-ends", "--run-metacluster",
        "--only-kmeans", "--input-reads", "--kmeans-number-of-rounds"};
    
    void change_parameter(string key, string value){
        if(find(names.begin(), names.end(), key) == names.end()){
            cerr << "Invalid parameter name: " << key << endl;
            exit(-1);
        }
        
        if(key == "--filtering-good-out"){
            good_reads_out = value; 
        } else if(key == "--filtering-bad-out"){
            bad_reads_out = value; 
        } else if(key == "--preclusters-out"){
            preclusters_out = value; 
        } else if(key == "--final-clusters-out"){
            final_clusters_out = value; 
        } else if(key == "--workspace"){
            workspace_path = value; 
        } else if(key == "--threads"){
            number_of_threads = stoi(value); 
        } else if(key == "--read-filtering-k"){
            read_filter_k = stoi(value); 
        } else if(key == "--read-filtering-threshold"){
            read_filter_threshold = stoi(value); 
        } else if(key == "--preclustering-k"){
            grouping_k = stoi(value); 
        } else if(key == "--precluster-min-size"){
            min_group_size = stoi(value); 
        } else if(key == "--precluster-max-size"){
            max_group_size = stoi(value); 
        } else if(key == "--kmeans-kmer-length"){
            k_means_k_mer_length = stoi(value); 
        } else if(key == "--kmeans-number-of-clusters"){
            k_means_number_of_clusters = stoi(value); 
        } else if(key == "--bitvector"){
            bitvector = value; 
        } else if(key == "--paired-ends"){
            paired_ends = read_bool(value); 
        } else if(key == "--only-kmeans"){
            only_kmeans = true;
        } else if(key == "--input-reads"){
            input_reads = value;
        } else if(key == "--run-metacluster"){
            run_metacluster = read_bool(value); 
        } else if(key == "--kmeans-number-of-rounds"){
            k_means_rounds = stoi(value);
        }  else{
            cerr << "Invalid parameter value: " << key << ": " << value << endl;
            exit(-1);
        }
    }
};





void print_config(Config config){
    cerr << "Config:" << endl;
    cerr << "===============" << endl;
    cerr << "Filtering (k, threshold): " << config.read_filter_k << ", " << config.read_filter_threshold << endl;
    cerr << "Preclustering (k, min group size, max group size): " << config.grouping_k << ", " << config.min_group_size << ", " << config.max_group_size << endl;
    cerr << "k-means: (k-mer length, number of clusters): " << config.k_means_k_mer_length << ", " << config.k_means_number_of_clusters << endl;
    cerr << "Number of threads: " << config.number_of_threads << endl;
    cerr << "Paired ends: " << (config.paired_ends ? "Yes" : "No") << endl;
    cerr << "Use Metacluster's k-means: " << (config.run_metacluster ? "Yes" : "No") << endl;
    cerr << "Bit vector type: " << config.bitvector << endl;
    cerr << "===============" << endl;
}



Config parse_config_file(string filename){
    
    ifstream in(filename);
    if(!in.good()){
        cerr << "Error: Could not open config file: " << filename << endl;
        exit(1);
    }
    
    Config C;
    set<string> parameters_read;
    while(true){
        string line; getline(in,line);
        if(in.eof()) break;
        trim(line);
        if(line == "" || line[0] == '#') continue; // '#' Means comment line
        stringstream ss(line);
        
        string key; ss >> key;
        
        // Trim whitespace and colon from key
        while(isspace(key.back()) || key.back() == ':'){
            key.pop_back();
            if(key.size() == 0){
                cerr << "Invalid config file (did you forget ':' after a parameter name?" << endl;
                exit(-1);
            }
        }
        
        string value; ss >> value;
        C.change_parameter(key,value);
        parameters_read.insert(key);
    }
    
    if(parameters_read.size() != C.names.size()){
        cerr << "Not all config parameters defined\nRequired parameters:\n";
        for(auto key : C.names){
            cerr << key;
            if(parameters_read.count(key) == 0){
                cerr << " <- MISSING";
            }
            cerr << endl;
        }
        exit(1);
    }
    return C;
}

void check_run_environment(Config config){
    
    // Check that the workspace is writable
    struct stat info;
    if(stat(config.workspace_path.c_str(), &info) != 0){
        cerr << "Error: Cannot access " << config.workspace_path << endl;
        cerr << strerror(errno) << endl;
        exit(1);
    } else if(!(info.st_mode & S_IFDIR)){
        cerr << "Error: " << config.workspace_path << " is not a directory" << endl;
        exit(1);
    }
    
    // Check that the input reads are readbale
    ifstream input(config.input_reads);
    if(!input.good()){
        cerr << "Error: cannot read file " << config.input_reads << endl;
        exit(1);
    }
    
    // Check that the output locations are writable
    ofstream goodOut, badOut, preclustersOut, finalclustersOut;
    goodOut.open(config.good_reads_out, ofstream::app);
    badOut.open(config.bad_reads_out, ofstream::app);
    preclustersOut.open(config.preclusters_out, ofstream::app);
    finalclustersOut.open(config.final_clusters_out, ofstream::app);
    
    if(!goodOut.good()){
        cerr << "Error: cannot write to " << config.good_reads_out << endl;
        exit(1);
    }
    if(!badOut.good()){
        cerr << "Error: cannot write to " << config.bad_reads_out << endl;
        exit(1);
    }
    if(!preclustersOut.good()){
        cerr << "Error: cannot write to " << config.preclusters_out << endl;
        exit(1);
    }
    if(!finalclustersOut.good()){
        cerr << "Error: cannot write to " << config.final_clusters_out << endl;
        exit(1);
    }
}

int64_t index_of(int64_t x, vector<int64_t>& v){
    // Binary search
    int64_t pos = 0;
    int64_t step = 1LL << 32;
    while(step != 0){
        if(pos + step < v.size() && v[pos+step] <= x)
            pos += step;
        step /= 2;
    }
    if(v[pos] == x) return pos;
    return -1; // Not found
}



// Update mate pair list to match the reads in the filtered set of reads
void update_mate_pairs(vector<int64_t>& new_ids, string paired_end_file){
    
    if(paired_end_file == "") return;
    sort(new_ids.begin(), new_ids.end());
    
    ifstream pairs_in(paired_end_file);
    ofstream pairs_out(paired_end_file + ".temp");
    while(true){
        string line;
        getline(pairs_in, line);
        if(!pairs_in.good()) break;
        stringstream ss(line);
        int64_t r1, r2;
        ss >> r1 >> r2;
        if(index_of(r1,new_ids) != -1 && index_of(r2, new_ids) != -1)
            pairs_out << index_of(r1,new_ids) << " " << index_of(r2, new_ids) << "\n";
    }
    copy_file(paired_end_file + ".temp", paired_end_file);
}



// run_all: NOTE: calls delete on C and kmeans
template<class t_bitvector> 
void compute_preclusters(BWT_cluster<t_bitvector>* C, Config config, string paired_end_file){
    
    cerr << getTimeString() << " Running file " << config.input_reads << endl;
    cerr << getTimeString() << " Starting" << endl;
    
    if(config.paired_ends){
        // Parse the data to generate the list of mate pairs... todo: do somewhere else
        ifstream input(config.input_reads);
        Parser parser;
        parser.parseConcatenate(input,'$', paired_end_file);
    }
    
    cerr << getTimeString() << " Removing duplicates" << endl;
    string no_duplicates = config.workspace_path + "/input_sorted_no_duplicates.txt";
    sort_remove_duplicates_update_mate_pairs(config.input_reads, no_duplicates, paired_end_file, paired_end_file + ".temp"); // does not update mate pairs if paired_end_file == ""
    if(config.paired_ends) copy_file(paired_end_file + ".temp", paired_end_file);
    
    C->loadData(no_duplicates);

    vector<int64_t> kept_reads = C->filterReads(config.read_filter_k,config.read_filter_threshold,true,config.good_reads_out, config.bad_reads_out);

    delete C;
    
    // Have to keep the read set sorted for ropebwt and the division of the BWT for parallelism
    sort_fasta(config.good_reads_out, config.good_reads_out + ".temp");
    copy_file(config.good_reads_out + ".temp", config.good_reads_out);
    if(config.paired_ends) update_mate_pairs(kept_reads, paired_end_file);
    
    BWT_cluster<sdsl::bit_vector>* C2 = new BWT_cluster<sdsl::bit_vector>(config.workspace_path, paired_end_file, config.number_of_threads);
    C2->loadData(config.good_reads_out);

    ofstream precluster_read_ids(config.workspace_path + "/precluster_read_ids.txt");
    C2->groupReads(config.grouping_k,true,config.max_group_size,precluster_read_ids);
    precluster_read_ids.close(); // flush

    ifstream precluster_read_ids_in(config.workspace_path + "/precluster_read_ids.txt");
    ofstream filtered_precluster_read_ids_out(config.workspace_path + "/precluster_read_ids_filtered.txt");;
    C2->filterGroups(precluster_read_ids_in,config.min_group_size,filtered_precluster_read_ids_out);
    filtered_precluster_read_ids_out.close(); // flush
    delete C2;
    
    sort_fasta(config.good_reads_out,config.good_reads_out); // To make the read ids match with the fasta file
    collect_clusters(config.good_reads_out,config.workspace_path + "/precluster_read_ids_filtered.txt", config.preclusters_out); // Map read ids into read strings


}

void run_kmeans(Config config){
    cerr << getTimeString() << " Starting k-means phase" << endl;
    if(config.run_metacluster){
        K_means_metacluster kmeans(config.preclusters_out, config.final_clusters_out, config.k_means_k_mer_length);
        kmeans.run(config.final_clusters_out,config.workspace_path);
    } else{
        K_means kmeans(config.preclusters_out, true, config.k_means_k_mer_length,config.number_of_threads);
        kmeans.run(config.k_means_number_of_clusters,config.k_means_rounds,config.final_clusters_out,config.workspace_path); // TODO: 100 to command line parameter
    }
}

int main(int argc, char** argv){

    
    
    if(argc % 2 == 0){
        cerr << "Usage: BWT_cluster --input-reads [filename.fna] (--config [your_config_file.txt])" << endl;
        cerr << "Config file defaults to ./config.txt" << endl;
        cerr << "The parameters of the pipeline are defined in the config file. The parameters can also be" << endl;
        cerr << "overridden by giving command line arguments. See config.txt for details" << endl;
        return -1;
    }
    
    std::time_t start_time_seconds = std::time(nullptr);

    string configFile = "./config.txt"; // The default config file
    
    // See if the user gave a custom configuration file
    for(int i = 1; i < argc-1; i++){
        if(string(argv[i]) == "--config"){
            configFile = string(argv[i+1]);
            break;
        }
    }
    Config config = parse_config_file(configFile);
    config.only_kmeans = false;
    
    // Read the command line parameters.
    // Parameters in given on the command line override the parameters in the config file
    for(int i = 1; i < argc-1; i++){
        if(string(argv[i]) == "--config"){
            i++; // Config file has already been parsed before
        } else{
            string key = argv[i];
            string value = argv[i+1];
            config.change_parameter(key,value);
            i++;
        }
    }
    
    print_config(config);
    check_run_environment(config);
    
    string paired_end_file = config.workspace_path + "/paired_reads.txt";
    if(config.paired_ends == false) paired_end_file = "";
    
    if(config.bitvector == "basic"){
        BWT_cluster<sdsl::bit_vector>* C = new BWT_cluster<sdsl::bit_vector>(config.workspace_path, paired_end_file, config.number_of_threads);
        if(!config.only_kmeans) compute_preclusters(C, config, paired_end_file);
        run_kmeans(config);
    }
    else if(config.bitvector == "rrr"){
        BWT_cluster<sdsl::rrr_vector<>>* C = new BWT_cluster<sdsl::rrr_vector<>>(config.workspace_path, paired_end_file, config.number_of_threads);      
        if(!config.only_kmeans) compute_preclusters(C, config, paired_end_file);
        run_kmeans(config);
    } else{
        cerr << "Error: Invalid bitvector. This should never happen" << endl;
        exit(1);
    }
    
    cerr << getTimeString() << " Peak memory: " << getPeakRSS() << " bytes" << endl;
    cerr << getTimeString() << " Total time: " << std::time(nullptr) - start_time_seconds << " seconds" << endl;
    cerr << getTimeString() << " Finished" << endl;
    
}

