/* 
    Read NCBI FASTA format files into an internally re-used buffer
    
    See http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
*/

#ifndef __fasta_reader_HPP
#define __fasta_reader_HPP
#include <istream>

class fasta_reader
{
private:
    /* pImpl pattern allows private members to be defined in cpp file */
    struct privates;
    privates * m;
public:
    /* Initialises the buffer and reads the header line
       Next() must be called before any of the Get methods
       Caller retains resonsibility for disposing of inputFile 
    */
    fasta_reader(std::istream * inputFile);
    
    /* Destructor */
    virtual ~fasta_reader();
    
    /* Read the next sequence into the given sequence buffer. Needs to be called to read the first sequence
       Returns false if we are at end of file. 
       Concatenates sequences that are spread across multiple lines, and strips out whitespace and numbers. 
       */
    bool next();
    
    /* Returns the current sequence, null-terminated 
       Returns null if Next() has not been called, or Next() returned false indicating end of file
       The pointer is valid until Next() is called again. */
    const char * get_sequence();
};

#endif
