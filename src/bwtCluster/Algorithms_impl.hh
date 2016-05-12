#include <sdsl/suffix_array_algorithm.hpp>
#include "Iterators.hh"
#include "BD_BWT_index.hh"
#include <thread>

template <class t_bitvector>
class Mark_right_maximal_kmers_algorithm{
private:
    BD_BWT_index<t_bitvector>* index;
    int64_t k;
    int64_t tau;
    int64_t max_number_of_threads;
    bool fill;
    atomic<int64_t> number_of_threads;
    std::mutex marking_mutex;
    std::mutex new_thread_mutex;
    sdsl::bit_vector* marks;
    
    void launch_thread(typename BD_BWT_index_iterator<t_bitvector>::Stack_frame frame){
        // Launch a new thread
        number_of_threads++;
        BD_BWT_index_iterator<t_bitvector> it_new(index);
        it_new.iteration_stack.clear();
        it_new.iteration_stack.push_back(frame);
        it_new.current = frame;
        thread t(&Mark_right_maximal_kmers_algorithm::iteration_thread, this, it_new);
        t.detach();
    }

    void iteration_thread(BD_BWT_index_iterator<t_bitvector> it){
        while(it.next(k)){
            new_thread_mutex.lock();
            if(it.iteration_stack.size() >= 10 && number_of_threads < max_number_of_threads){
                auto frame = it.iteration_stack.front();
                it.iteration_stack.pop_front();
                launch_thread(frame);
            }
            new_thread_mutex.unlock();
            
            Interval& I = it.current.intervals.forward;
            if(I.size()  >= tau){
                // Fill interval with ones
                marking_mutex.lock();
                if(fill){
                    for(int64_t i = I.left; i <= I.right; i++)
                        (*marks)[i] = 1;
                    
                } else{
                    (*marks)[I.left] = 1;
                    (*marks)[I.right] = 1;
                }
                marking_mutex.unlock();
            }
        }
        
        number_of_threads--;
    }
    
    void launch_thread_rc(typename BD_BWT_index_RC_iterator<t_bitvector>::Stack_frame frame){
        // Launch a new thread
        number_of_threads++;
        BD_BWT_index_RC_iterator<t_bitvector> it_new(index);
        it_new.iteration_stack.clear();
        it_new.iteration_stack.push_back(frame);
        it_new.current = frame;
        thread t(&Mark_right_maximal_kmers_algorithm::iteration_thread_rc, this, it_new);
        t.detach();
    }

    void iteration_thread_rc(BD_BWT_index_RC_iterator<t_bitvector> it){
        while(it.next(k)){
            new_thread_mutex.lock();
            if(it.iteration_stack.size() >= 10 && number_of_threads < max_number_of_threads){
                auto frame = it.iteration_stack.front();
                it.iteration_stack.pop_front();
                launch_thread_rc(frame);
            }
            new_thread_mutex.unlock();
            
            Interval& I_normal = it.current.intervals.forward;
            Interval& I_rc = it.current.intervals_rc.forward;
            int64_t size;
            if(I_normal.left == I_rc.left){
                size = I_normal.size(); // k-mer is the rc of itself
            } else{
                size = I_normal.size() + I_rc.size();
            }
            
            if(size >= tau){ 
                marking_mutex.lock();
                
                // Mark normal interval
                for(int64_t i = I_normal.left; i <= I_normal.right; i++)
                    (*marks)[i] = 1;
                // Mark rc interval
                for(int64_t i = I_rc.left; i <= I_rc.right; i++)
                    (*marks)[i] = 1;
                
                marking_mutex.unlock();
            }
        }
        
        number_of_threads--;
    }
    
public:
    Mark_right_maximal_kmers_algorithm(BD_BWT_index<t_bitvector>& index, int64_t k, int64_t tau, bool fill, int64_t nThreads) :
        index(&index), k(k), tau(tau), max_number_of_threads(nThreads), fill(fill), number_of_threads(0) {}
        
    sdsl::bit_vector run_non_rc(){
        sdsl::bit_vector bitvec(this->index->forward_bwt.size(),0); // Initialize to zeroes
        marks =  &bitvec;
        BD_BWT_index_iterator<t_bitvector> it(this->index);
        number_of_threads = 1;
        thread t(&Mark_right_maximal_kmers_algorithm::iteration_thread, this, it);
        t.detach();
        while(number_of_threads > 0) std::this_thread::sleep_for(std::chrono::milliseconds(10));
        return bitvec;
    }
    
    sdsl::bit_vector run_rc(){
        sdsl::bit_vector bitvec(this->index->forward_bwt.size(),0); // Initialize to zeroes
        marks =  &bitvec;
        BD_BWT_index_RC_iterator<t_bitvector> it(this->index);
        number_of_threads = 1;
        thread t(&Mark_right_maximal_kmers_algorithm::iteration_thread_rc, this, it);
        t.detach();
        while(number_of_threads > 0) std::this_thread::sleep_for(std::chrono::milliseconds(10));
        return bitvec;
    }
};

template <class t_bitvector>
sdsl::bit_vector get_right_maximal_kmers(BD_BWT_index<t_bitvector>& index, int64_t k, int64_t tau, bool fill, int64_t nThreads){
    Mark_right_maximal_kmers_algorithm<t_bitvector> algorithm(index,k,tau,fill,nThreads);
    return algorithm.run_non_rc();
}

template <class t_bitvector>
sdsl::bit_vector get_right_maximal_kmers_rc(BD_BWT_index<t_bitvector>& index, int64_t k, int64_t tau, int64_t nThreads){
    // Here the fill parameter is unused, so I arbitrarily set it to true
    Mark_right_maximal_kmers_algorithm<t_bitvector> algorithm(index,k,tau,true,nThreads);
    return algorithm.run_rc();
}

namespace find_kmer_rc_pairs_namespace{
    
    // Common data shared between all threads
    template <class t_bitvector>
    class Common_data{
    public:
        atomic<int64_t> number_of_threads;
        int64_t max_number_of_threads;
        int64_t k;
        std::mutex new_thread_mutex;
        std::mutex marking_mutex;
        sdsl::bit_vector* marks;
        std::ostream* out;
        BD_BWT_index<t_bitvector>* index;
    };
    
    template <class t_bitvector>
    void iteration_thread(BD_BWT_index_RC_iterator<t_bitvector> it, Common_data<t_bitvector>& common){
        while(it.next(common.k)){
            common.new_thread_mutex.lock();
            if(it.iteration_stack.size() >= 10 && common.number_of_threads < common.max_number_of_threads){
                // Launch a new thread
                common.number_of_threads++;
                auto frame = it.iteration_stack.front();
                it.iteration_stack.pop_front();
                BD_BWT_index_RC_iterator<t_bitvector> it_new(common.index,frame);
                thread t(&iteration_thread<t_bitvector>, it_new, ref(common));
                t.detach();
            }
            common.new_thread_mutex.unlock();
            
            // Each position will be marked twice, because both it and its RC interval
            // will be iterated. However it makes no difference.
            Interval_pair& I = it.current.intervals;
            Interval_pair& I_rc = it.current.intervals_rc;
            if(I == I_rc) continue; // kmer is a reverse-complement of itself
            if(I.forward.size() >= 1 && I_rc.reverse.size() >= 1){
                common.marking_mutex.lock();
                (*common.marks)[I.forward.left] = 1; // An arbitrary position from the interval
                (*common.marks)[I_rc.forward.left] = 1; // An arbitrary position from the interval
                *common.out << I.forward.left << " " << I_rc.forward.left << "\n";
                common.marking_mutex.unlock();
            }
        }
        common.number_of_threads--;
    }
};

template <class t_bitvector>
sdsl::bit_vector find_kmer_rc_pairs(BD_BWT_index<t_bitvector>& index, int64_t k, int64_t nThreads, std::ostream& out){
    find_kmer_rc_pairs_namespace::Common_data<t_bitvector> common;
    sdsl::bit_vector marks(index.forward_bwt.size(),0); // Initialize to zeroes
    common.marks = &marks;
    common.out = &out;
    common.index = &index;
    common.number_of_threads = 1;
    common.max_number_of_threads = nThreads;
    common.k = k;
    
    BD_BWT_index_RC_iterator<t_bitvector> root(&index);
    thread t(&find_kmer_rc_pairs_namespace::iteration_thread<t_bitvector>, root, ref(common));
    t.detach();
    while(common.number_of_threads > 0) std::this_thread::sleep_for(std::chrono::milliseconds(10));
    return marks;
}

void fill_intervals_length_at_least_tau(sdsl::bit_vector& v, int64_t tau){

    // Fill the intervals that have length at least tau with ones and zero out of rest
    int64_t interval_start = -1;
    int64_t state = 0; // 0 = not inside an interval, 1 = inside an interval
    for(int64_t i = 0; i < v.size(); i++){
        if(state == 0 && v[i] == 0); // Outside an interval: do nothing
        else if(state == 0 && v[i] == 1){ // Entering an interval
            state = 1;
            interval_start = i;
        }
        else if(state == 1 && v[i] == 0) v[i] = 1; // Inside an interval
        else if(state == 1 && v[i] == 1){ // Closing an interval
            if(i - interval_start + 1 < tau){ // Oops, too short, have to zero it out
                for(int64_t j = interval_start; j <= i; j++){
                    v[j] = 0;
                }
            }
            state = 0;
        }
    }
}


