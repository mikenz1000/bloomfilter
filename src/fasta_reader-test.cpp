#include "unit_test.hpp"
#include "fasta_reader.hpp"
#include <fstream>
#include <string>

class test_reader_t : public unit_test
{
    void operator() ()
    {
        section("reading fasta file");
        std::ifstream f("data/test_long.fa");
        fasta_reader r(&f);
        check(r.next(), "fasta_reader: first Next() call");
        std::string desired1("GAACGGTCCGGCCGCATCCATTTCTTCCCTGTAGCGAATCGCGAAAATCGTCCGGAGTCTTAGTGTCTAAAGGTGGTTCACACGGAGATATGAGCGCGCC");
        std::string desired2("GCTAGGGCTTCCGGTTCCATGTGGAGATTAGCCCGTAAATTCAGATCCCGATCGCCACCACTCGATCAACTTTTCTTTAACATGTATTCATTTCCAATAA");
        check(std::string(r.get_sequence()) == desired1, "fasta_reader: first sequence read");
        check(r.next(), "FASTReader: second Next() call");
        check(std::string(r.get_sequence()) == desired2, "fasta_reader: second sequence read");
    }
    
} test_reader;