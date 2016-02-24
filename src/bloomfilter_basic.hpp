/*
    A basic bloom filter implementation that stores the bits in blocks of a given type (block_t)
    
    block_t must be unsigned.  Recommended values are uint64_t, uint32_t depending on architecture
    
    index_t is the type used as an index for the set/test operations
    
    The hash function can be provided, or otherwise defaults to the standard std::hash

    See BloomFilter.hpp for explanation of the methods
*/
#ifndef __bloomfilter_basic_HPP
#define __bloomfilter_basic_HPP
#include "bloomfilter.hpp"

template<typename index_t, typename block_t, typename Hash = std::hash<index_t> >
class bloomfilter_basic : public bloomfilter<index_t>
{
protected:
    block_t * bitarray;   
    size_t blockcount;
public:     

    bloomfilter_basic(int m, int h) : bloomfilter<index_t>(m,h)
    {
        // round up so that if we for example have m=9 and sizeof(block_t)=8 then we get 2 elements in the array
        // (9+1*8-1)/(1*8) = 16/8 = 2
        blockcount = (m+sizeof(block_t)*8-1)/(sizeof(block_t)*8);
        bitarray = new block_t[blockcount];
        clear();
    };

    virtual void set(const index_t & kmer)
    {
        const size_t BitsPerElement = sizeof(block_t) * 8;
        static Hash hashfunction;
        index_t hashvalue = kmer;
        for (int hcount = this->h; hcount > 0; hcount--)
        {
            hashvalue = hashfunction(hashvalue);
            // we expect the compiler to automatically turn this into a shift because it's a const power of two
            size_t bitindex = hashvalue % this->m;
            size_t offset = bitindex / BitsPerElement;
            block_t mask = ( ((block_t)1) << (bitindex & (BitsPerElement-1)));
            bitarray[offset] |= mask;
        }
    }

    virtual bool test(const index_t & kmer) const
    {
        const size_t BitsPerElement = sizeof(block_t) * 8;
        static Hash hashfunction;
        index_t hashvalue = kmer;
        for (int hcount = this->h; hcount > 0; hcount--)
        {
            hashvalue = hashfunction(hashvalue);
            size_t bitindex = hashvalue % this->m;
            // we expect the compiler to automatically turn this into a shift because it's a const power of two
            size_t offset = (bitindex) / BitsPerElement;
            block_t mask = ( ((block_t)1) << (bitindex & (BitsPerElement-1)));
            if (!(bitarray[offset] & mask)) return false;
        }
        return true;
    }
    
    void clear()
    {
        memset(bitarray, 0, blockcount*sizeof(block_t));
    }

    virtual ~bloomfilter_basic()
    {
        delete bitarray;
    }
};

#endif