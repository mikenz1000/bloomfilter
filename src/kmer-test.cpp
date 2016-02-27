#include "unit_test.hpp"
#include "kmer.hpp"
#include "fasta_reader.hpp"
#include <sstream>

class test_kmer_t : public unit_test
{
    void operator() ()
    {
        section("reading kmer from fasta file");
        std::stringstream f(">SEQ1\nGGATA\n");
        fasta_reader r(&f);
        check(r.next(), "first Next() call");
        kmer_ops ops(3);
        const char * i = r.get_sequence();
        kmer_t kmer;
        check(ops.read_first(&kmer,&i), "read first kmer");
        check(ops.str(kmer) == "GGA", "first kmer is correct");
        check(ops.read_next(&kmer,&i), "read 2nd kmer");
        check(ops.str(kmer) == "GAT", "2nd kmer is correct");
        check(ops.read_next(&kmer,&i), "read 3rd kmer");
        check(ops.str(kmer) == "ATA", "3rd kmer is correct");
        check(!ops.read_next(&kmer,&i), "end read");
    }
} test_kmer;
