#ifndef ITERATORS_HH
#define ITERATORS_HH

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include "BD_BWT_index.hh"
#include <algorithm>

template<typename t_bitvector> class BD_BWT_index;


/**
 * Class BWT_inversion_iterator
 * 
 * An iterator that iterates backwards over a continuous range in the text from which the
 * given BWT has been build.
 * 
 */
template<class t_bitvector>
class BWT_inversion_iterator{
private: 
    
    // char_block_starts[c] = the lexicographical rank of the first lexicographical
    // suffix to start with character c
    vector<int64_t> char_block_starts; // For backward step
    
    vector<uint8_t> alphabet;
    int64_t done;
    
public:
    sdsl::wt_huff<t_bitvector>* bwt;
    int64_t pos; // current bwt position
    int64_t last; // last bwt position to be iterated, inclusive
    uint8_t current_char;
    
    // The read id of the read that contains the suffix at pos. In
    // case the suffix starts with a dollar, read_id is the id of the read that
    // follows the dollar
    int64_t read_id; 
   
    /** @brief Go to the next text position
     * 
     *  Updates the member variables pos, current_char and read_id
     * 
     *  @return false iff the current position is the last position, i.e. there is not "next" node
     */
    bool next();

     /**
     * 
     * Construct a new iterator. NOTE: suffixes "from" and "to" MUST start with a dollar
     * 
     *  @param from      bwt position of the suffix FOLLOWING the first suffix to be iterated
     *  @param to        bwt position of last suffix to be iterated (inclusive)
     *  @param bwt       a wavelet tree of the bwt to be inverted
     *  @param read_id   the read id of the first read to be iterated
     * 
     */
    BWT_inversion_iterator(sdsl::wt_huff<t_bitvector>* bwt, int64_t from, int64_t to, int64_t read_id) :
        char_block_starts(255), done(false), bwt(bwt),
        pos(from), last(to), current_char('$'), read_id(read_id+1) {
            alphabet = get_interval_symbols(*bwt, Interval(0,bwt->size()-1));
            sort(alphabet.begin(), alphabet.end());
            count_smaller_chars(*bwt,alphabet,char_block_starts,Interval(0,bwt->size()-1));
        }
    
};

/**
 * Class BD_BWT_index_iterator
 * 
 * Iterates the suffix link tree of the given index.
 * 
 */
template<class t_bitvector>
class BD_BWT_index_iterator{
    
public:
    
    class Stack_frame{
    public:
        Interval_pair intervals;  // forward interval, reverse interval
        int64_t depth; // depth in the suffix link tree
        char extension; // the label on the arc between this node and its parent
        Stack_frame(Interval_pair intervals, int64_t depth, char extension) : intervals(intervals), depth(depth), extension(extension) {}
        Stack_frame(){}
    };
    
    BD_BWT_index<t_bitvector>* index;
    
    // Iteration state
    std::deque<Stack_frame> iteration_stack;
    Stack_frame current; // A stack frame containing information about the current node
    std::string label; // The string on the path from the root to the current node
    
    // Reused space between iterations
    std::vector<int64_t> char_counts;
    
    BD_BWT_index_iterator(BD_BWT_index<t_bitvector>* index);
    
    // Initialize to the given stack frame. The label will be starting from the depth of f
    BD_BWT_index_iterator(BD_BWT_index<t_bitvector>* index, Stack_frame f);
    
    /**
     * @brief Go to the next suffix link tree node
     * 
     * Updates the member variables iteration_stack, current and label.
     * 
     * @return False if all nodes have been iterated, else true
     */
    bool next(); // Go to the next node. Return false if none found
    
    /**
     * @brief Go to the next depth k suffix link tree node
     * 
     * Updates the member variables iteration_stack, current and label.
     * 
     * @return False if all nodes have been iterated, else true
     */
    bool next(int64_t k);
    
private:
    void push_right_maximal_children(Stack_frame f);
    void update_label(Stack_frame f);
};

/**
 * Class BD_BWT_index_RC_iterator
 * 
 * Iterates the generalized suffix link tree of the string corresponding to the
 * given index, concatenated with its reverse complement
 * 
 */
template<class t_bitvector>
class BD_BWT_index_RC_iterator{
public:
    
    class Stack_frame{
    public:
        Interval_pair intervals;  // forward interval, reverse interval
        Interval_pair intervals_rc;   // forward interval, reverse interval 
        int64_t depth; // depth in the suffix link tree
        char extension; // the label on the arc between this node and its parent in the forward string
        Stack_frame(Interval_pair intervals, Interval_pair intervals_rc, int64_t depth, char extension) 
            : intervals(intervals), intervals_rc(intervals_rc), depth(depth), extension(extension) {}
        Stack_frame(){}
    };
    
    BD_BWT_index<t_bitvector>* index;
    
    // Iteration state
    std::deque<Stack_frame> iteration_stack;
    Stack_frame current; // A stack frame containing information about the current node
    std::string label; // The string on the path from the root to the current node
    
    // Reused space between iterations:
    std::vector<int64_t> char_counts;
    std::vector<int64_t> char_counts_rc;
    std::vector<uint8_t> child_right_extensions;
    std::vector<uint8_t> child_rc_left_extensions;
    std::vector<uint64_t> ranks_i_unused; // For the get_interval_symbols call
    std::vector<uint64_t> ranks_j_unused; // For the get_interval_symbols call
    std::vector<bool> GST_extensions; // Generalized suffix tree (string + RC of string) extensions

    // Initialize to the root node
    BD_BWT_index_RC_iterator(BD_BWT_index<t_bitvector>* index);
    // Initialize to the given stack frame. The label will be starting from the depth of f
    BD_BWT_index_RC_iterator(BD_BWT_index<t_bitvector>* index, Stack_frame f);

    
    /**
     * @brief Go to the next suffix link tree node
     * 
     * Updates the member variables iteration_stack, current and label.
     * 
     * @return False if all nodes have been iterated, else true
     */
    bool next(); // Go to the next node. Return false if none found
    
    /**
     * @brief Go to the next depth k suffix link tree node
     * 
     * Updates the member variables iteration_stack, current and label.
     * 
     * @return False if all nodes have been iterated, else true
     */
    bool next(int64_t k);
    
private:
    void push_right_maximal_children(Stack_frame f);
    void update_label(Stack_frame f);
};


#include "Iterators_impl.hh"

#endif