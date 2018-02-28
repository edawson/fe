#ifndef H_FE_INC
#define H_FE_INC
#include <cstdint>
#include <sstream>
#include <string.h>
#include "tabix.hpp"
#include "tbx.h"
#include "regidx.h"
#include "pliib.hpp"


struct vcf_compact{
    string line = "";
    string chrom;
    int64_t pos;
    int64_t end;
    vector<string> other_cols;
    vector<string> info_keys;
    vector<string> info_vals;
};


// Parse a VCF info section to a pair of vectors
pair<vector<string>, vector<string>> info_to_vecs(string info);

// convert a VCF line (string) to a simple struct
vcf_compact parse_vcf_line(char* line);

// Read a VCF file kinda fast


// Chunk a VCF
// return an array of chunks
vector<vector<char*> > chunk_a_vcf(string filename);

// annoatate a single variant if it has a tabix annotation
//void annotate(vcf_compact& vc);




//string annotate(string& chrom, int& pos, int& end, string& infos, Tabix t);
string annotate(string& chrom, int& pos, int& end, string& infos, tbx_t* tt, htsFile* fp);
//string annotate(string& line, tbx_t* tt, htsFile* fp);

// Annotate a chunk of variants
//void annotate(vector<string>& vlines, vector<vcf_compact>& annos);

#endif
