#include <vector>
#include "BD_BWT_index.hh"
#include "Iterators.hh"
#include <iostream>

using namespace std;

template<class t_bitvector>
void BD_BWT_index_iterator<t_bitvector>::push_right_maximal_children(Stack_frame f){
    count_smaller_chars(index->forward_bwt,index->alphabet,char_counts,f.intervals.forward);
    for(uint8_t c : index->alphabet){
        if(c == '$') continue;
        Interval_pair child = left_extend(*index,f.intervals,c,char_counts[c]);
        if(child.forward.size() == 0) continue; // Extension not possible
        if(is_right_maximal(*index,child)){
            // Add child to stack
            iteration_stack.push_back(Stack_frame(child,f.depth+1,c));
        }
    }    
}

template<class t_bitvector>
void BD_BWT_index_iterator<t_bitvector>::update_label(Stack_frame f){
    while(label.size() > 0 && label.size() >= f.depth) // Unwind stack
        label.pop_back();
    if(f.extension != 0)
        label.push_back(current.extension);    
}

template<class t_bitvector>
bool BD_BWT_index_iterator<t_bitvector>::next(int64_t k){
    
    while(true){
        if(iteration_stack.empty()) return false;
        
        current = iteration_stack.back(); iteration_stack.pop_back();
        update_label(current);
        
        if(current.depth == k) return true; // Stop recursing to children and give control back

        push_right_maximal_children(current);

    }
}

template<class t_bitvector>
bool BD_BWT_index_iterator<t_bitvector>::next(){
    if(iteration_stack.empty()) return false;
    
    current = iteration_stack.back(); iteration_stack.pop_back();
    update_label(current);
    push_right_maximal_children(current);
    
    return true;
}


inline char complement(char c){
    if(c == 'A') return 'T';
    if(c == 'C') return 'G';
    if(c == 'G') return 'C';
    if(c == 'T') return 'A';
    if(c == '$') return '$';
    assert(false); // Should not come here
}

// In the generalized suffix tree:
// Let S be the set of right-extensions of the forward string and
// S_{RC} be the set of complements of left-extensions of the RC string. A string is
// right-maximal iff |S union S_{RC}| >= 2
template<class t_bitvector>
void BD_BWT_index_RC_iterator<t_bitvector>::push_right_maximal_children(Stack_frame f){
    count_smaller_chars(index->forward_bwt,index->alphabet,char_counts,f.intervals.forward);
    count_smaller_chars(index->reverse_bwt,index->alphabet,char_counts_rc,f.intervals_rc.reverse);
    for(char c : index->alphabet){
        
        if(c == '$' || c == 0) continue;
        Interval_pair child = left_extend(*index,f.intervals,c,char_counts[c]);
        if(child.forward.size() == 0) continue;
        Interval_pair child_rc = right_extend(*index,f.intervals_rc,complement(c),char_counts_rc[complement(c)]);
        
        sdsl::int_vector_size_type nExtensions, nExtensions_rc;
        get_interval_symbols<t_bitvector>(index->reverse_bwt, child.reverse, nExtensions, 
                                            child_right_extensions, ranks_i_unused, ranks_j_unused);
        get_interval_symbols<t_bitvector>(index->forward_bwt, child_rc.forward, nExtensions_rc,
                                            child_rc_left_extensions, ranks_i_unused, ranks_j_unused);
        for(auto c : index->alphabet) GST_extensions[c] = false; // Clear
        for(int64_t i = 0; i < nExtensions; i++) GST_extensions[child_right_extensions[i]] = true;
        for(int64_t i = 0; i < nExtensions_rc; i++) GST_extensions[complement(child_rc_left_extensions[i])] = true;
        int64_t nDistinct = 0;
        for(auto c : index->alphabet) if(GST_extensions[c]) nDistinct++;
        
        // The intervals in child can not be empty, but the intervals in child_rc might be empty
        if(nDistinct >= 2 || ((child.reverse.size() + child_rc.forward.size() >= 2) && GST_extensions['$'] == true)){
            // Two distinct extensions or two or more dollars
            // Add child to stack to be iterated on later
            iteration_stack.push_back(Stack_frame(child, child_rc, f.depth+1, c));
        }
    }
}

template<class t_bitvector>
void BD_BWT_index_RC_iterator<t_bitvector>::update_label(Stack_frame f){
    while(label.size() > 0 && label.size() >= f.depth) // Unwind stack
        label.pop_back();
    if(f.extension != 0)
        label.push_back(current.extension);    
}

template<class t_bitvector>
bool BD_BWT_index_RC_iterator<t_bitvector>::next(){
    if(iteration_stack.empty()) return false;
    
    current = iteration_stack.back(); iteration_stack.pop_back();
    update_label(current);       
    push_right_maximal_children(current);
    
    return true;
}


template<class t_bitvector>
bool BD_BWT_index_RC_iterator<t_bitvector>::next(int64_t k){
    while(true){
        if(iteration_stack.empty()) return false;
        
        current = iteration_stack.back(); iteration_stack.pop_back();
        update_label(current);
        
        if(current.depth == k) return true; // Stop recursing to children and give control back
        
        push_right_maximal_children(current);

    }
}

template<class t_bitvector>
bool BWT_inversion_iterator<t_bitvector>::next(){
    if(done) {
        return false;   
    }
    
    if(current_char == '$') read_id--;
    
    // Do a backward step  
    char prev_char = (*bwt)[pos];
    int64_t new_pos = char_block_starts[prev_char] + bwt->rank(pos, prev_char); 
    if(new_pos == last) done = true; // This iteration is ok, but all the following will return false
    pos = new_pos;
    current_char = prev_char;
    
    return true;
}



template<class t_bitvector>
BD_BWT_index_iterator<t_bitvector>::BD_BWT_index_iterator(BD_BWT_index<t_bitvector>* index) : index(index), char_counts(256) {
    Interval empty_string(0,index->forward_bwt.size()-1);
    iteration_stack.push_back(Stack_frame(Interval_pair(empty_string,empty_string), 0, 0));
    current = iteration_stack.back();
    label = "";
}

template<class t_bitvector>
BD_BWT_index_iterator<t_bitvector>::BD_BWT_index_iterator(BD_BWT_index<t_bitvector>* index, Stack_frame f) : index(index), char_counts(256) {
    Interval empty_string(0,index->forward_bwt.size()-1);
    iteration_stack.push_back(f);
    current = f;
}

template<class t_bitvector>
BD_BWT_index_RC_iterator<t_bitvector>::BD_BWT_index_RC_iterator(BD_BWT_index<t_bitvector>* index) 
: index(index), char_counts(256), char_counts_rc(256), 
  child_right_extensions(index->forward_bwt.sigma), child_rc_left_extensions(index->forward_bwt.sigma), 
  ranks_i_unused(index->forward_bwt.sigma), ranks_j_unused(index->forward_bwt.sigma), 
  GST_extensions(256) {
    Interval empty_string = Interval(0,index->forward_bwt.size()-1);
    iteration_stack.push_back(Stack_frame(Interval_pair(empty_string,empty_string),Interval_pair(empty_string,empty_string), 0, 0));
    current = iteration_stack.back();
    label = "";
}

template<class t_bitvector>
BD_BWT_index_RC_iterator<t_bitvector>::BD_BWT_index_RC_iterator(BD_BWT_index<t_bitvector>* index, Stack_frame f) : 
  index(index), char_counts(256), char_counts_rc(256), 
  child_right_extensions(index->forward_bwt.sigma), child_rc_left_extensions(index->forward_bwt.sigma), 
  ranks_i_unused(index->forward_bwt.sigma), ranks_j_unused(index->forward_bwt.sigma), 
  GST_extensions(256){
    
    iteration_stack.push_back(f);
    current = f;
}
