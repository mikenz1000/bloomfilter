/*
    A portable bloom filter implementation that just uses vector<bool> to store the bit array
    
    In fact, this can perform well depending on the platform..
    
    See BloomFilter.hpp for explanation of the arguments/methods
*/
#ifndef __BLOOMFILTER_BITSET_HPP
#define __BLOOMFILTER_BITSET_HPP
#include "bloomfilter.hpp"
#include <vector>

template<typename index_t, typename Hash = std::hash<index_t> >
class bloomfilter_vectorbool : public bloomfilter<index_t>
{
protected:
    /* m bits of storage */
    std::vector<bool> bitarray;
public:     

    bloomfilter_vectorbool(int m, int h) : bloomfilter<index_t>(m,h), bitarray(m,false)
    {
    };

    void set(const index_t & kmer)
    {
        Hash hashfunction = Hash();
        index_t hashvalue = kmer;
        for (int hcount = this->h; hcount > 0; hcount--)
        {
            hashvalue = hashfunction(hashvalue);
            size_t offset = hashvalue % this->m;
            bitarray[offset] = true;
        }
    }

    bool test(const index_t & kmer) const
    {
        Hash hashfunction = Hash();
        index_t hashvalue = kmer;
        for (int hcount = this->h; hcount > 0; hcount--)
        {
            hashvalue = hashfunction(hashvalue);
            if (!bitarray[hashvalue % this->m]) return false;
        }
        return true;
    }
    
    void clear()
    {
        bitarray.assign(this->m,false);
    }

};

#endif