/*
    A cheating bloom filter that uses an unordered set so never returns a false positive
    
    See BloomFilter.hpp for explanation of the methods
*/
#ifndef __BLOOMFILTER_PERFECTCHEAT_HPP
#define __BLOOMFILTER_PERFECTCHEAT_HPP
#include "bloomfilter.hpp"
#include <unordered_set>

template<typename index_t, typename Hash = std::hash<index_t> >
class bloomfilter_perfectcheat : public bloomfilter<index_t>
{
protected:
    std::unordered_set<index_t> storage;
public:     

    bloomfilter_perfectcheat(int m, int h) : bloomfilter<index_t>(m,h)
    {
    };

    void set(const index_t & kmer)
    {
        storage.insert(kmer);
    }

    bool test(const index_t & kmer) const
    {
        return storage.count(kmer) > 0;
    }

    void clear()
    {
        storage.clear();
    }

};

#endif