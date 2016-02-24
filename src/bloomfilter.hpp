/* 
  Bloom Filter abstract base class (a few different implementations will be compared)
  
  The bloom filter acts as a memory efficient set that may return false positives.
  
  If the hash argument is not supplied you get the default C++ hash
  
  index_t is the type of the index parameter that will be used to do the set/test operations
  
  Algorithm is explained in theh introductory parts of this paper:
  Space-efficient and exact de Bruijn graph representation based on a Bloom filter
      http://almob.biomedcentral.com/articles/10.1186/1748-7188-8-22
*/
#ifndef __BLOOMFILTER_HPP
#define __BLOOMFILTER_HPP
#include <functional>
#include <cmath>

template<typename index_t>
class bloomfilter
{
public: 
    /* Reset the bloom filter - you can assume it is cleared when it is first created */
    virtual void clear() = 0;
    
    /* Add the given kmer to the set */
    virtual void set(const index_t & kmer) = 0;
    
    /* Return true if we think kmer is in the set (may return false positives) */
    virtual bool test(const index_t & kmer) const = 0;
    
    /* Destructor */
    virtual ~bloomfilter() {};   
    
    /* Given n the number of items inserted, returns the expected false positive probability 
        (1 - exp(-kn/m))^k    
    */
    double expected_false_positive_probability(double n)
    {
        return pow(1-exp(-h * n / m), h);
    }
    
    /* Given p the desired false probability and n the expected number of elements, determine the required size 
       assuming optimum h is used */
    static size_t determine_m(double p, double n)
    {
        return - n * log(p) / (log(2)*log(2));
    }
    
    /* Given m the filter size and n the expected number of elements, return optimum number of hash functions */
    static size_t determine_h(double m, double n)
    {
        return std::ceil(m / n * log(2));
    }
    
    int getm() { return m; }
    int geth() { return h; }
protected:
    /* Constructor - m is the filter array size in bits, h is the number of hash functions */
    bloomfilter(int m, int h = 0) : m(m), h(h) {};
    
    /* the number of bits */
    int m;
    
    /* the number of hash functions */
    int h;
};

#endif