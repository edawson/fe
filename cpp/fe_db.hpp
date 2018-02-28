#ifndef H_FE_DB
#define H_FE_DB
#include <string>
#include <sstream>
#include "pliib.hpp"
#include <omp.h>
#include "tbx.h"
#include "kstring.h"
#include "htsfile.h"

using namespace std;
// Defines a class for managing tab-delimited text files
// to be used as input to the annotation pipeline
class FE_DB{
    private:
        tbx_t * tabix_db;

    public:
        FE_DB(string filename);
        ~FE_DB();
        string get_header();
        //int number_of_variants();
    
        void query(const char* region, kstring_t& ret);
        //void query(string region, kstr& ret);
        //void query(char* chrom, int pos, int end, kstr& ret);
};
#endif
