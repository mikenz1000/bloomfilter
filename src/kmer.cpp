#include <stdexcept>
#include "kmer.hpp"

struct text_mapping_t {
    /* Stores the unshifted NUCLEOTIDE_ values for the given ascii characters or kmer_invalid otherwise */
    kmer_t asciiToBits[128];
    char bitsToAscii[4];
    text_mapping_t()
    {
        // generate the mapping from ascii to bits
        for (int i = 0;i < 128;i ++) asciiToBits[i] = KMER_INVALID;
        asciiToBits[(size_t)'A'] = NUCLEOTIDE_A;
        asciiToBits[(size_t)'C'] = NUCLEOTIDE_C;
        asciiToBits[(size_t)'G'] = NUCLEOTIDE_G;
        asciiToBits[(size_t)'T'] = NUCLEOTIDE_T;
        
        // and the reverse
        bitsToAscii[NUCLEOTIDE_A] = 'A';
        bitsToAscii[NUCLEOTIDE_C] = 'C';
        bitsToAscii[NUCLEOTIDE_G] = 'G';
        bitsToAscii[NUCLEOTIDE_T] = 'T';        
    }
} text_mapping;

kmer_ops::kmer_ops(kmer_size_t length) : length(length)
{
    if (length > (sizeof(kmer_t)*4-1)) throw std::runtime_error("kmer length too long");
    
    // build the mask
    mask = 0;
    for (int i = 0;i < length;i ++)
    {
        mask <<= 2;
        mask |= 3;
    }
}

bool kmer_ops::read_first(kmer_t * kmer, char const * * sequence) const
{
    kmer_size_t todo = length;
    *kmer = 0;
    while (todo-- > 0) 
    {
        // load next nucleotide
        if (!read_next(kmer,sequence)) return false;  
    }
    return true;
}

bool kmer_ops::read_next(kmer_t * kmer,  char const * * sequence) const
{
    // read the next char
    const char * i = *sequence;   
    size_t ix = (size_t)*i;
    if (!ix) return false;
    if (ix > 127) throw std::runtime_error("Non-ascii characters encountered when reading nucleotides");
    (*sequence)++;
    
    // map to bits
    kmer_t bits = text_mapping.asciiToBits[ix];
    if (bits == KMER_INVALID) throw std::runtime_error("Invalid character when reading nucleotides");
    
    // shift into working copy
    *kmer <<= 2;
    *kmer |= bits;
    
    // superfluous when called from ReadFirst but leaving as-is
    *kmer &= mask;
    return true;
}

std::string kmer_ops::str(const kmer_t & kmer) const
{
    // allocate the string
    std::string buf(length, '?');
    
    // and write from end to start
    kmer_t working = kmer;
    kmer_size_t todo = length;
    while (todo-- > 0)
    {
        buf[todo] = text_mapping.bitsToAscii[working & 3];
        working >>= 2;
    }
    return buf;
}

kmer_t kmer_ops::complement(const kmer_t kmer) const
{
    return kmer ^ mask;
}