#include "fe.hpp"
#include "tabix.hpp"
#include <iostream>

using namespace std;
int main(int argc, char** argv){

    // Read in VCF
    string st  = string (argv[1]);
    Tabix qt(st);

    string dbstr = string(argv[2]);
    tbx_t* tt = tbx_index_load(dbstr.c_str());
    htsFile* fp = hts_open(dbstr.c_str(), "r");

    string empt = "";
    string one = "1";
    int pos = 1;
    int end = 10000000;
    // Open DB tabix iterator
    // For each variant (one day  in parallel )
    string line;
    while (qt.getNextLine(line)){
        vector<string> splits = pliib::split(line, '\t');
        int pos = stoi(splits[1]);
        int end = stoi(splits[1]);
        annotate(splits[0], pos, end, splits[7], tt, fp);
    }
    // Annotate it with annotations from each annotation file
    // just pass the info column to the annotator function
    // Output that vcf

    return 0;
}
