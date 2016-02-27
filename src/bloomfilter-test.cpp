#include <iostream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <unordered_set>
#include <algorithm>
#include "unit_test.hpp"
#include "terminal.hpp"
#include "scoped_timer.hpp"
#include "bloomfilter.hpp"
#include "bloomfilter_basic.hpp"
#include "bloomfilter_vectorbool.hpp"
#include "bloomfilter_perfectcheat.hpp"
#ifdef USESSE
#include "bloomfilter_sse.hpp"
#endif
#include "kmer.hpp"

// quick test of basic functionality
class test_quick_t : public unit_test 
{
    public:
    void operator()()
    {
        section("quick bloom filter test");
        bloomfilter_basic<uint64_t,kmer_t> bf(100,5);
        bf.set(1);
        bf.set(2);
        check(bf.test(1), "TestBloomFilter: Test 1");
        check(bf.test(2), "TestBloomFilter: Test 2");
        bf.set(64);
        bf.set(65);
        check(bf.test(64), "TestBloomFilter: Test 1");
        check(bf.test(65), "TestBloomFilter: Test 2");
    }  
} test_quick;

// in-depth test
class test_speed_t : public unit_test 
{
public:
    // the number of insertions we will make
    int n_count = 200000;
    
    // how many times to we repeat the fill and test operations in order
    // to get a decent timing
    int repeat = 100;
    
    // the size of the filter in bits
    int m;

    // the number of hash functions
    int h;

    // intensive test and timing
    template<typename T>
    void test_bloomfilter(const char * info)
    {
        // output a description of the filter
        std::cout << info << ": ";
        
        // create the bloom filter
        T bf(m,h);
        
        // build our data set
        std::minstd_rand0 rng (243345);  // minstd_rand0 is a standard linear_congruential_engine
        std::unordered_set<kmer_t> data;
        for (int i = 0;i < n_count;i ++)
            data.insert((kmer_t) rng());
        
        // time the fill operation
        {
            scoped_timer t("\tfilling", n_count);
            for (int i = 0;i < repeat;i ++)
                std::for_each(data.begin(), data.end(), [&bf](const kmer_t &i){bf.set(i);});
        }

        // time the test operation for all the values we know we set
        {
            scoped_timer t("\ttesting", n_count);
            int false_negative = 0;
            for (int i = 0;i < repeat;i ++)
            {
                std::for_each(data.begin(), data.end(), [&bf,&false_negative](const kmer_t &i)
                {
                    if (!bf.test(i)) false_negative++;
                });
            }
            if (false_negative) 
                std::cout << terminal::red << "\nFalse negative rate: " << (false_negative / (double)n_count) << terminal::reset << std::endl;
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
        std::cout << "\tp(false +ve) " << std::fixed << (false_positive / (double)(done)) << 
            " vs E[p(false +ve)] " << bf.expected_false_positive_probability(done) << terminal::reset << std::endl;
        
    }

    void operator() ()
    {        
        double p = 0.01;
        while (p <= 0.01)
        {
            m = bloomfilter<kmer_t>::determine_m(p,n_count);
            h = bloomfilter<kmer_t>::determine_h(m,n_count);
            std::cout << "Determined m = " << m << ", h = " << h << " for desired p(false +ve) " << p << std::endl;

        #ifdef USESSE
            bloomfilter_sse<kmer_t> bfsse_cpphash(m,h);
            test_bloomfilter(&bfsse_cpphash, "sse C++ hash");
        #endif   
        
        #ifdef TIME_KMER_HASH
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,0> >("64 bit blocks                     ");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,1> >("64 bit blocks - 1 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,2> >("64 bit blocks - 2 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,3> >("64 bit blocks - 3 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,4> >("64 bit blocks - 4 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,4> >("64 bit blocks - 5 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,4> >("64 bit blocks - 6 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,4> >("64 bit blocks - 7 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,kmer_hash,4> >("64 bit blocks - 8 byte mis-aligned");

           test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t,kmer_hash,0> >("32 bit blocks                     ");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t,kmer_hash,1> >("32 bit blocks - 1 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t,kmer_hash,2> >("32 bit blocks - 2 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t,kmer_hash,3> >("32 bit blocks - 3 byte mis-aligned");
            
           test_bloomfilter< bloomfilter_basic<kmer_t,uint8_t,kmer_hash> >   ("8 bit blocks                      ");
        #endif   
        
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,0> >("64 bit blocks, std::hash...       ");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,1> >("64 bit blocks - 1 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,2> >("64 bit blocks - 2 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,3> >("64 bit blocks - 3 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,4> >("64 bit blocks - 4 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,4> >("64 bit blocks - 5 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,4> >("64 bit blocks - 6 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,4> >("64 bit blocks - 7 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint64_t,std::hash<kmer_t>,4> >("64 bit blocks - 8 byte mis-aligned");

           //test_bloomfilter< bloomfilter_vectorbool<kmer_t, kmer_hash> >             ("vector<bool>, kmer hash           ");           
           test_bloomfilter< bloomfilter_vectorbool<kmer_t> >                        ("vector<bool>, std::hash           ");     
           
           test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t,std::hash<kmer_t>,0> >("32 bit blocks, std::hash...       ");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t,std::hash<kmer_t>,1> >("32 bit blocks - 1 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t,std::hash<kmer_t>,2> >("32 bit blocks - 2 byte mis-aligned");
           test_bloomfilter< bloomfilter_basic<kmer_t,uint32_t,std::hash<kmer_t>,3> >("32 bit blocks - 3 byte mis-aligned");
            
           test_bloomfilter< bloomfilter_basic<kmer_t,uint8_t,std::hash<kmer_t> > >  ("8 bit blocks                      ");

            p*=10;
        }
    }
} test_timing;