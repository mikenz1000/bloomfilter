/*
   kmer is a value type that stores a kmer of a length known only to be less than MAXKMERLENGTH
   
   kmer_ops is a host class instantiated with the length parameter.  It provides the methods for
   manipulating the kmers.  This way we don't have to store the kmer length in the kmer value type
   or pass it explicitly every time.  Also allows kmer_ops to cache handy bitmasks and things
*/
#ifndef __KMER_HPP
#define __KMER_HPP

#include <string>
#include "../thirdparty/spookyhash/SpookyV2.h"

/* the kmer typedef */
#ifndef MAXKMERLENGTH
#error MAXKMERLENGTH must be defined (suggest 31)
#endif
#if MAXKMERLENGTH < 32
// we don't use the full 32 because then we can't distinguish between KMER_INVALID and a potentially valid kmer
typedef uint64_t kmer_t;
#else
#error MAXKMERLENGTH too long (can't be more than 64)
#endif

/* the unsigned type to be used for length and index values */
typedef size_t kmer_size_t;

/* the invalid value for a kmer - the use of the top 2 bits makes it impossible for a valid kmer to have the same value */
const kmer_t KMER_INVALID = (kmer_t)-1;

/* the best hash function we know of for kmers */
struct kmer_hash
{
    kmer_hash() {}
    kmer_hash(const kmer_hash & src) {}
    std::size_t operator()(kmer_t const& kmer) const 
    {
         return (size_t)SpookyHash::Hash64(&kmer, sizeof(kmer_t), 0);
    }
};

/* the way the dna letters are encoded 
   chosen so that the complement is achieved by reversing the bits */
#define NUCLEOTIDE_A 0
#define NUCLEOTIDE_C 1
#define NUCLEOTIDE_G 2
#define NUCLEOTIDE_T 3

/* the kmewOps class 
   I see a need for a kmer length variable that I don't want to pass
   around with every kmer.  Rather than establish a global variable this kmer_ops class records the length
   and hosts the methods that manipulate kmers
*/
class kmer_ops
{
private:
    kmer_size_t length;  
    
    /* used when shifting left to mask out bits that have gone out the top 
       also used to generate complements */
    kmer_t mask;
public:
    kmer_ops(kmer_size_t length);
    
    /* advances sequence to the next character after the first kmer
       returns false if the sequence was not long enough to generate a single kmer */
    bool read_first(kmer_t * kmer, char const * * sequence) const;
    
    /* updates kmer and advances sequence . returns false if we ran out and we can't generate a new kmer */
    bool read_next(kmer_t * kmer,  char const * * sequence) const;
    
    /* renders the kmer to text */
    std::string str(const kmer_t & kmer) const;
    
    /* generates the complement A->T, C->G, G->C, T->A */
    kmer_t complement(const kmer_t kmer) const;
}; 

#endif