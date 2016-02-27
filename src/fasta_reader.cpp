#include "fasta_reader.hpp"
#include <vector>
#include <algorithm>
#include <stdexcept>


#include <iostream>

/* The maximum anticipated line size including terminating null - if exceeded it causes exceptions */
const int MAXLINESIZE = 140;
    
/* Private members and functions */
struct fasta_reader::privates
{
    std::vector<char> header;
    std::vector<char> buffer;
    std::istream * inputFile;
};

/* Public */
fasta_reader::fasta_reader(std::istream * inputFile)
{
    m = new privates();
    m->inputFile = inputFile;
    m->header.resize(MAXLINESIZE);
    m->buffer.resize(MAXLINESIZE);
}
    
fasta_reader::~fasta_reader()
{
    delete m;
}
    
bool fasta_reader::next()
{
    // ensure the next character is >
    char gt = m->inputFile->get();
    
    // check for end of file
    if (gt == EOF) return false; 
    
    // validate we have the beginning of sequence marker
    if (gt != '>') throw std::runtime_error("Expected > character but read something else");
    
    // read until the end of the line
    m->inputFile->getline(&m->header[0], m->header.size());
    if (m->inputFile->fail()) throw std::runtime_error("Sequence header longer than maximum line size");
    
    // now we can start reading the sequence
    size_t writePos = 0;
    m->inputFile->getline(&m->buffer[writePos], m->buffer.size());
    if (m->inputFile->fail()) throw std::runtime_error("Buffer line longer than maximum allowed line size");
    
    // peak to see if we have more
    gt = m->inputFile->peek();
    while ( (gt != '>') && (gt != EOF) )
    {
        // move writePos to the occurance of the null (that must exist as per getLine documentation)
        writePos = std::find(m->buffer.begin()+writePos, m->buffer.end(), 0) - m->buffer.begin();
        
        // expand the buffer
        m->buffer.resize(writePos + MAXLINESIZE);
        
        // read
        m->inputFile->getline(&m->buffer[writePos], MAXLINESIZE);
        if (m->inputFile->fail()) throw std::runtime_error("Subsequent buffer line longer than maximum allowed line size");

        // and peek for the > on the next line
        gt = m->inputFile->peek();
    }
    return true;
}
    
const char * fasta_reader::get_sequence()
{
    return &m->buffer[0];
}