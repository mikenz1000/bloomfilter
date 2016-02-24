/*
    An experiment coding the bloom filter using SSE instructions (128 bit blocks)
    
    See BloomFilter.hpp for explanation of the arguments/methods
*/
#ifndef __BLOOMFILTER_SSE_HPP
#define __BLOOMFILTER_SSE_HPP
#include "bloomfilter.hpp"
#include <xmmintrin.h>
#include <stdlib.h>
#include <iostream>

template<typename index_t, typename Hash = std::hash<index_t> >
class bloomfilter_sse : public bloomfilter<index_t>
{
protected:
    typedef __m128 block_t;
    /* m bits of storage */
    block_t * bitarray;
    size_t blockcount;
    uint8_t * storage;
    __m128 * masks;
public:     

    bloomfilter_sse(int m, int h) : bloomfilter<index_t>(m,h)
    {
        const size_t BitsPerElement = sizeof(block_t) * 8;

        // round up so that if we for example have m=9 and sizeof(block_t)=8 then we get 2 elements in the array
        // (9+1*8-1)/8 = 16/8 = 2
        size_t bits = (m+sizeof(block_t)*8-1);
        blockcount = bits / 128;
        
        // we also at the same time allocate the space for the precalculated bit masks
        size_t bytes = 128 * sizeof(block_t) + blockcount * sizeof(block_t);
        storage = (uint8_t*)malloc(bytes);
        
        // advance the pointer in order to get 16 byte alignment (couldn't find aligned malloc etc on this mac)
        masks = (__m128*)(storage + (((size_t)storage) & 15)); 
        bitarray = masks + 128;

        // precalculate the bit masks
        __m128 __attribute__ ((aligned (16))) one = _mm_set_epi64x(1, 0);
        for (int i = 0;i < 128;i ++)
        {
            masks[i] = one;
            one = _mm_slli_si128(one, 1);
        }
    };

    void set(const index_t & kmer)
    {
        const size_t BitsPerElement = sizeof(block_t) * 8;
        const Hash hashfunction = Hash();
        kmer_t hashvalue = kmer;
        for (int hcount = this->h; hcount > 0; hcount--)
        {
            hashvalue = hashfunction(hashvalue);
            // we expect the compiler to automatically turn this into a shift because it's a const power of two
            size_t offset = (hashvalue % this->m) / BitsPerElement;
            size_t maskindex = hashvalue & (BitsPerElement-1);
            bitarray[offset] = _mm_or_ps(bitarray[offset], masks[maskindex]);
        }
    }

    bool test(const index_t & kmer) const
    {
        __m128 __attribute__ ((aligned (16))) zero = _mm_setzero_si128();
        const size_t BitsPerElement = sizeof(block_t) * 8;
        const Hash hashfunction = Hash();
        kmer_t hashvalue = kmer;
        for (int hcount = this->h; hcount > 0; hcount--)
        {
            hashvalue = hashfunction(hashvalue);
            // we expect the compiler to automatically turn this into a shift because it's a const power of two
            size_t offset = (hashvalue % this->m) / BitsPerElement;
            if (_mm_movemask_epi8(
                _mm_cmpeq_epi32(
                _mm_and_ps(bitarray[offset],masks[hashvalue & (BitsPerElement-1)]),zero)) != 0xFFFF) return false;
        }
        return true;
    }

    void clear()
    {
        memset(bitarray, 0, this->blockcount * sizeof(block_t));
    }

    virtual ~bloomfilter_sse()
    {
        std::cout << "freeing ... ";
        free(storage);
        std::cout << "... done" << std::endl;
    }
};

#endif