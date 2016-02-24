/*
    This is a primitive unit test archicture - I have developed more elegant ones in my professional work!
*/
#include <iostream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <unordered_set>
#include <algorithm>

#include "../thirdparty/spookyhash/SpookyV2.h"
#include "scoped_timer.hpp"
#include "bloomfilter.hpp"
#include "bloomfilter_basic.hpp"
#include "bloomfilter_vectorbool.hpp"
#include "bloomfilter_perfectcheat.hpp"
#ifdef USESSE
#include "bloomfilter_sse.hpp"
#endif

// colour codes for terminal output
const std::string red("\033[0;31m");
const std::string green("\033[1;32m");
const std::string yellow("\033[1;33m");
const std::string cyan("\033[0;36m");
const std::string magenta("\033[0;35m");
const std::string reset("\033[0m");

// the type for the bloom filter index
typedef uint64_t kmer_t;

/* another hash to try */
struct kmer_hash
{
    kmer_hash() {}
    kmer_hash(const kmer_hash & src) {}
    std::size_t operator()(kmer_t const& kmer) const 
    {
         return (size_t)SpookyHash::Hash64(&kmer, sizeof(kmer_t), 0);
    }
};

// simple assertion function writes to stderr
void test(bool result, const std::string & message)
{
    if (!result) 
        std::cerr << "FAILED: " << message << std::endl;
}

// quick test of basic functionality
void test_quick()
{
    bloomfilter_basic<uint64_t,kmer_t> bf(100,5);
    bf.set(1);
    bf.set(2);
    test(bf.test(1), "TestBloomFilter: Test 1");
    test(bf.test(2), "TestBloomFilter: Test 2");
    bf.set(64);
    bf.set(65);
    test(bf.test(64), "TestBloomFilter: Test 1");
    test(bf.test(65), "TestBloomFilter: Test 2");
}

// the number of insertions we will make
int n_count = 1000000;

// the size of the filter in bits
int m;

// the number of hash functions
int h;

// intensive test and timing
template<typename T>
void test_bloomfilter(const char * info)//bloomfilter<kmer_t> * bf, const char * info)
{
    // output a description of the filter
    std::cout << info << ":\n";
    
    // create the bloom filter
    T bf(m,h);
    
    // build our data set
    std::minstd_rand0 rng (243345);  // minstd_rand0 is a standard linear_congruential_engine
    std::unordered_set<kmer_t> data;
    for (int i = 0;i < n_count;i ++)
        data.insert((kmer_t) rng());
    
    // time the fill operation
    {
        scoped_timer t("\tFilling", n_count);
        std::for_each(data.begin(), data.end(), [&bf](const kmer_t &i){bf.set(i);});
    }

    // time the test operation for all the values we know we set
    {
        scoped_timer t("\tTesting", n_count);
        int false_negative = 0;
        std::for_each(data.begin(), data.end(), [&bf,&false_negative](const kmer_t &i)
        {
            if (!bf.test(i)) false_negative++;
        });
        if (false_negative) 
            std::cout << red << "\tFalse negative rate: " << (false_negative / (double)n_count) << reset << std::endl;
    }
        
    // Test for false positives by additional random numbers ...
    int false_positive = 0;
    int todo = n_count;
    int done = 0;
    for (int i = 0;i < todo; i++)
    {
        kmer_t k = rng();
        bool in_bf = bf.test(k);
        bool in_data = data.count(k) > 0;
        // ... and ignoring them if they were already in the set
        if (!in_data) 
        {
            done ++;
            if (in_bf) false_positive ++;
        } 
    }
    std::cout << green << "\tFalse positive rate " << std::fixed << (false_positive / (double)(done)) << 
        " versus expected " << bf.expected_false_positive_probability(done) << reset << std::endl;
      
}

int main(int argc, const char* args[])
{
    test_quick();
    
    double p = 0.001;
    while (p <= 0.1)
    {
        m = bloomfilter<kmer_t>::determine_m(p,n_count);
        h = bloomfilter<kmer_t>::determine_h(m,n_count);
        std::cout << "Determined m = " << m << ", h = " << h << " for desired p(false +ve) " << p << std::endl;

    #ifdef USESSE
        bloomfilter_sse<kmer_t> bfsse_cpphash(m,h);
        test_bloomfilter(&bfsse_cpphash, "sse C++ hash");
    #endif   
//        test_bloomfilter< bloomfilter_perfectcheat<kmer_t> >("Perfect cheat");
//        test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t> >("64 bit blocks");
//        test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t> >("32 bit blocks");
//        test_bloomfilter< bloomfilter_basic<kmer_t,uint16_t> >("16 bit blocks");
//        test_bloomfilter< bloomfilter_basic<kmer_t,uint8_t> >("8 bit blocks");
        test_bloomfilter< bloomfilter_vectorbool<kmer_t> >("vector<bool>");     
        test_bloomfilter< bloomfilter_vectorbool<kmer_t, kmer_hash> >("vector<bool>, kmer hash");
        
        p*=10;
    }
    return 0;
}