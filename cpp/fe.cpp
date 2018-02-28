#include "fe.hpp"

using namespace std;
// Parse a VCF info section to a pair of vectors
pair<vector<string>, vector<string>> info_to_vecs(string info){

}

void get_tbx_header(string& header){
    
}

// convert a VCF line (string) to a simple struct
vcf_compact parse_vcf_line(string line){
    vcf_compact vc;
    vector<string> tokens = pliib::split(line, '\t');
    vc.line = line;
    //vector<string> i_tokens = split(tokens[7], ';');
    //for (int i = 0; i < i_tokens.size(); i++){
    //      
    //}
}



// Chunk a VCF
// return a vector of header lines and an array of chunks
/**
 * tuple<vector<string>, vector<vector<string> > > chunk_a_vcf(string filename){
    vector<string> header;
    vector<vector<string> > lines
    vector<string> chunk;
    ifstream infile;
    infile.open(filename);
    string line;
    int count = 0;
    int limit = 1000;
    while ( getline(infile, line) ){
        if (line[0] == '#'){
            header.push_back(line);
        }
        else{
           chunk.push_back(line);
           count += 1;
           if (count == limit){
                lines.push_back(chunk);
                chunk.clear();
                count = 0;
           }
        }
    }
    return make_tuple(header, lines);
}
*/

//string annotate(string& chrom, int& pos, int& end, string& infos, Tabix tabi){
string annotate(string& chrom, int& pos, int& end, string& infos, tbx_t* tt, htsFile* fp){
    vector<string> splits = pliib::split(infos, ';');
    stringstream st;
    st << chrom << ":" << pos << "-" << end;
    const char* reg = st.str().c_str();
    
    //tbx_itr_queryi(tt, 1, pos, end);
    hts_itr_t* iter = tbx_itr_querys(tt, reg);
    kstring_t str = {0,0,0};
    int nseq;
    const char** seq;
    seq = tbx_seqnames(tt, &nseq);
    if (iter){
        while(tbx_itr_next(fp, tt, iter, &str) >= 0){
            cout << str.s << endl << str.l << endl;
                //cout << reg  << "=" << seq[iter->curr_tid] << ", " <<
                //pos << "=" << iter->curr_beg << ", " << end << "=" << iter->curr_end << endl << infos << endl;
            if (!(strcmp(reg, seq[iter->curr_tid]) == 0 || pos >= iter->curr_beg || end < iter->curr_end)) break;
        }
    }
    return "";
}

// Annotate a chunk of variants
//void annotate(vector<string>& vlines, vector<vcf_compact>& annos){
//
//}
