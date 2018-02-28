#include "fe_db.hpp"

using namespace std;

FE_DB::FE_DB(string filename){        
    tbx_t * tabix_db = tbx_index_load(filename.c_str());
    htsfile* fp = hts_open(filename.c_str(), "r");
}

FE_DB::~FE_DB(){
    delete tabix_db;
    delete fp;
}

string FE_DB::get_header(){
    stringstream st;
    kstring_t ks;
    while ( hts_getline(fn, KS_SEP_LINE, &ks) >= 0 ) {
        if ( !str.l || str.s[0]!=tbx->conf.meta_char ) {
            break;
        } else {
            st << string(ks.s) << "\n";
        }
    }
    return st.str();
}


void FE_DB::query(const char* region, char* ret, int& retlen){
    #pragma omp atomic
    {
        hts_itr_t* iter = tbx_itr_querys(tabix_db, region);
    }
    kstring_t str = {0,0,0};
    int nseq;
    const char** seq;
    seq = tbx_seqnames(tt, &nseq);
    if (iter){
        while(tbx_itr_next(fp, tt, iter, &str) >= 0){
            if (str.l > 0){
                ret = str.s;
                retlen = str.l;
            }
            if (!(strcmp(reg, seq[iter->curr_tid]) == 0 || pos >= iter->curr_beg || end < iter->curr_end)) break;
        }
    }

}

