#include "Read_group_builder.hh"
#include "Algorithms.hh"

template <class t_bitvector>
void Read_group_builder<t_bitvector>::write_groups(UnionFind& UF, std::ostream &output){
    vector<pair<int64_t,int64_t> > groups;
    for(int64_t r = 0; r < nReads; r++){
        int64_t r_cluster = UF.find(r);
        groups.push_back(make_pair(r,r_cluster));
    }

    // Sort by cluster id
    std::sort(groups.begin(), groups.end(), [](pair<int64_t,int64_t> a, pair<int64_t,int64_t> b){return a.second < b.second;});

    // Scan through the sorted list and write to output
    for(int64_t i = 0; i < groups.size(); i++){
        if(i != 0 && (groups[i].second != groups[i-1].second)) output << "\n";
        output << groups[i].first << " ";
    }

    output << "\n";
}

template<class t_bitvector>
void Read_group_builder<t_bitvector>::compute_marked_position_read_ids(BWT_inversion_iterator<t_bitvector> iter, sdsl::bit_vector& marks,
            sdsl::rank_support_v<1>& marks_rs, sdsl::int_vector<0>& read_ids, std::mutex& mutex){
     // Invert the BWT and collect the read ids
    while(iter.next()){
        if(marks[iter.pos] == 1){
            mutex.lock();
            read_ids[marks_rs.rank(iter.pos)] = iter.read_id;
            mutex.unlock();
        }
    }
}

template <class t_bitvector>
sdsl::bit_vector Read_group_builder<t_bitvector>::compute_k_submaximal_intervals(int64_t k){
    sdsl::bit_vector marks(index->forward_bwt.size(),0); // Initialize to zeroes
    
    BD_BWT_index_iterator<t_bitvector> it_k(index);
    while(it_k.next(k)){    
        Interval& I = it_k.current.intervals.forward;
        marks[I.left] = 1;
        marks[I.right] = 1;
    }
    
    sdsl::rank_support_v<1> marks_rs(&marks);
    sdsl::bit_vector marks_new = marks;
        
    BD_BWT_index_iterator<t_bitvector> it(index);
    Interval ancestor(-1,-1);
    while(it.next()){
        auto sf = it.current;
        if(sf.depth == k) {
            if(get_interval_symbols(index->forward_bwt, sf.intervals.forward).size() > 1)
                ancestor = Interval(-1,-1);
            else
                ancestor = sf.intervals.forward;
        }
        if(sf.depth > k){
            int64_t left = sf.intervals.forward.left;
            if(ancestor != Interval(-1,-1) && marks_rs.rank(left+1) % 2 == 1){
                // Current interval is a subinterval of another marked interval
                marks_new[ancestor.left] = 0;
                marks_new[ancestor.right] = 0;
            }
            if(get_interval_symbols(index->forward_bwt, sf.intervals.forward).size() > 1){
                // Left-maximal
                ancestor = Interval(-1,-1);
            }
        }
    }
    
    // Print statistics
    int64_t nRightMaximal = 0; int64_t nSubmaximal = 0;
    for(int64_t i = 0; i < index->forward_bwt.size(); i++){
        nRightMaximal += marks[i];
        nSubmaximal += marks_new[i];
    }
    cerr << getTimeString() << " Right-maximal " << k << "-mers: " << nRightMaximal << ", " << k << "-submaximal repeats: " << nSubmaximal << endl;
    return marks_new;
}

template <class t_bitvector>
UnionFind Read_group_builder<t_bitvector>::parallel_non_RC_union_find(int64_t maxGroupSize, sdsl::bit_vector& intervals){
    sdsl::rank_support_v<1> intervals_rs(&intervals);
    int64_t nIntervals = intervals_rs.rank(index->forward_bwt.size()) / 2;
    sdsl::int_vector<0> interval_reps(nIntervals,nReads,(int)log2(nReads) + 1);
    cerr << getTimeString() <<  " Interval representatives array takes " << size_in_bytes(interval_reps)
        << " bytes (" << size_in_bytes(interval_reps)*8 << " bits)" << endl;
        
    UnionFind UF; UF.init(nReads);
    
    
    // Join paired end reads first
    if(paired_end_file != ""){
        cerr << getTimeString() <<  " Joining paired end reads" << endl;
        ifstream pairs(paired_end_file);
        string line;
        while(getline(pairs,line)){
            stringstream ss(line);
            int64_t r1, r2;
            ss >> r1 >> r2;
            UF.doUnion(r1,r2);
        }
    }
    
    vector<thread> threads(nThreads);
    for(int64_t t = 0; t < nThreads; t++){
        thread x(&Read_group_builder::do_unions, this, iterators[t], ref(interval_reps), ref(UF), ref(intervals), ref(intervals_rs), maxGroupSize);
        threads[t] = move(x);
    }
    
    for(auto& x : threads) x.join();
    return UF;
}

template <class t_bitvector>
void Read_group_builder<t_bitvector>::do_unions(BWT_inversion_iterator<t_bitvector> iter,
                                                sdsl::int_vector<0>& interval_reps, UnionFind& UF, sdsl::bit_vector& intervals, 
                                                sdsl::rank_support_v<1>& intervals_rs, int64_t maxGroupSize){
    
    int64_t fragsize = 100; // From metacluster. Merge two preclusters if either is of size less than this
    while(iter.next()){
        int64_t rank = intervals_rs.rank(iter.pos + 1);
        if(rank % 2 == 1 || (rank % 2 == 0 && intervals[iter.pos] == 1)){
            // Inside an interval
            int64_t block = (rank - 1) / 2;
            
            interval_reps_mutex.lock();
            if(interval_reps[block] == nReads){// First time visiting this block
                interval_reps[block] = iter.read_id;
                interval_reps_mutex.unlock();
                continue;
            }
            interval_reps_mutex.unlock();
            
            // Note, the two lines below are a bit thread-unsafe, since the root of the
            // group containing a read might change in the middle of the find-operation.
            // However, it is not a problem for this algorithm, since if two reads are in
            // the same group originally, they will still be in the same group after
            // the changes.
            int64_t id1 = UF.find(interval_reps[block]);
            int64_t id2 = UF.find(iter.read_id);
            
            if(id1 != id2){
                union_find_mutex.lock();
                            
                int64_t size1 = UF.getSize(id1);
                int64_t size2 = UF.getSize(id2);
                if((size1 <= maxGroupSize && size2 <= maxGroupSize) || size1 < fragsize || size2 < fragsize)
                    UF.doUnion(id1,id2);
                if(size2 < size1)
                    interval_reps[block] = id2;

                union_find_mutex.unlock();
            }
        }
    }    
}

template <class t_bitvector>
void Read_group_builder<t_bitvector>::groupReads(int64_t k, int64_t maxGroupSize, std::ostream &output){
    cerr << getTimeString() << " Grouping reads" << endl;
    cerr << getTimeString() << " Computing the right-maximal interval bitvector (" << nThreads << " threads)" << endl;
    
    sdsl::bit_vector intervals = get_right_maximal_kmers(*index, k, 2, false, nThreads);
    
    cerr << getTimeString() << " Inverting the BWT and issuing union commands" << endl;
    UnionFind UF = parallel_non_RC_union_find(maxGroupSize,intervals);
    write_groups(UF,output); 
    
}

template <class t_bitvector>
void Read_group_builder<t_bitvector>::groupReadsSubmaximal(int64_t k, int64_t maxGroupSize, std::ostream &output){
    cerr << getTimeString() << " Grouping reads (k-submaximal)" << endl;    
    cerr << getTimeString() << " Computing the k-submaximal interval bitvector  NOT PARALLELLIZED" << endl;
    // Find the maximal intervals
    sdsl::bit_vector intervals = compute_k_submaximal_intervals(k);
    
    cerr << getTimeString() << " Inverting the BWT and issuing union commands" << endl; 

    UnionFind UF = parallel_non_RC_union_find(maxGroupSize,intervals);
    write_groups(UF,output);  
}

template <class t_bitvector>
void Read_group_builder<t_bitvector>::groupReadsRC(int64_t k, int64_t maxGroupSize, std::string workspace_path, std::ostream &output){
    cerr << getTimeString() << " Grouping reads (considering reverse complements)" << endl;
    cerr << getTimeString() << " Computing the right-maximal interval bitvector (" << nThreads << " threads)" << endl;
    
    // Mark the right-maximal k-mers and do the regular grouping
    sdsl::bit_vector intervals = get_right_maximal_kmers(*index, k, 2, false, nThreads);
    
    cerr << getTimeString() << " Joining the groups using non-RC matches" << endl;
    UnionFind UF = parallel_non_RC_union_find(maxGroupSize,intervals);
    
    // Find the right-maximal RC-kmers (kmers considering a string and its reverse complement the same string)
    // For each right-maximal RC-kmer, if it appears in both the forward and the RC string, mark one position
    // from the forward interval and one position from the reverse interval, and write the position pair to disk.
    // When done with all such kmers, reserve an integer array of length equal to the number of marks. Then, invert
    // the forward BWT, and when the current position is marked, write the read ID of the position
    // to the aforementioned interger array into the position determined by the rank of the mark. 
    // Then, stream the merge position pairs from disk, and fetch the read ID's of the pairs by rank 
    // queries to the integer array, and issue the Union-Find commands.
    
    cerr << getTimeString() << " Marking forward-RC pairs to be joined" << endl;;
    ofstream position_pairs_out(workspace_path + "/position_pairs.txt");
    sdsl::bit_vector marks = find_kmer_rc_pairs<t_bitvector>(*index,k,nThreads,position_pairs_out);
    sdsl::rank_support_v<1> marks_rs(&marks);
    int64_t nMarks = marks_rs.rank(marks.size());
    
    sdsl::int_vector<0> read_ids(nMarks,0,(int)log2(nReads) + 1);  
    
    cerr << getTimeString() << " Inverting the BWT to find the read ids of the marked positions" << endl;
    std::mutex read_id_mutex;
    vector<thread> threads(nThreads);
    for(int64_t t = 0; t < nThreads; t++){
        thread x(&Read_group_builder::compute_marked_position_read_ids, this, iterators[t], ref(marks), ref(marks_rs), ref(read_ids), ref(read_id_mutex));
        threads[t] = move(x);
    }
    
    for(auto& x : threads) x.join();
    
    cerr << getTimeString() << " Streaming the position pairs from disk and issuing Union commands" << endl;
    ifstream position_pairs_in(workspace_path + "/position_pairs.txt");
    assert(nMarks % 2 == 0);
    for(int64_t i = 0; i < nMarks/2; i++){
        int64_t pos1, pos2;
        position_pairs_in >> pos1 >> pos2;
        int64_t read1 = read_ids[marks_rs.rank(pos1)];
        int64_t read2 = read_ids[marks_rs.rank(pos2)];
        if(UF.getSize(read1) + UF.getSize(read2) <= maxGroupSize)
            UF.doUnion(read1, read2);
    }
    
    write_groups(UF,output); 
}