
#include "unit_test.hpp"

// // the type for the bloom filter index
// typedef uint64_t kmer_t;

// /* another hash to try */
// struct kmer_hash
// {
//     kmer_hash() {}
//     kmer_hash(const kmer_hash & src) {}
//     std::size_t operator()(kmer_t const& kmer) const 
//     {
//          return (size_t)SpookyHash::Hash64(&kmer, sizeof(kmer_t), 0);
//     }
// };

int main(int argc, const char* args[])
{
    unit_test::run_all();
}

