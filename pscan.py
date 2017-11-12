import sys
import multiprocessing
import subprocess
import tabix
import pyensembl
from cStringIO import StringIO
import argparse


bedtools = "bedtools"
bcftools = "bcftools"
bgzip = "bgzip"

tabixdb = None
bcsq_path = None
samtools_path = None

default_tb_database = ""
default_cpg_list = ""

tb_db = None

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest = "infile", required = True, type = argparse.FileType('r'),
        help = "A VCF file to sort, filter, slice and annotate")
parser.add_argument("-b", "--bed", dest = "bedfile", required = False, type = argparse.FileType('r'),
        help = "A BED file to restrict variant annotation to.", default = None)
parser.add_argument("-l", "--list", dest = "gene_list", required = False, type = argparse.FileType('r'),
        help = "A file containing a list of genes, one per line, to restrict variant annotation to.", default = None)
parser.add_argument("-t", "--tabix", dest = "tabix", required = False, type = argparse.FileType('r'),
        help = "A tabix-indexed and bgzip'ed database of variant annotations", default = None)
parser.add_argument("-c", "--consequence", dest = "csq", required = False, action = "store_true", help = "Run bcftools csq")
args = parser.parse_args()

header = ["#chr","pos(1-based)","ref","alt","aaref","aaalt",
"rs_dbSNP150","hg19_chr","hg19_pos(1-based)",
"hg18_chr","hg18_pos(1-based)","genename",
"cds_strand","refcodon","codonpos","codon_degeneracy",
"Ancestral_allele","AltaiNeandertal","Denisova","Ensembl_geneid",
"Ensembl_transcriptid","Ensembl_proteinid","aapos","SIFT_score",
"SIFT_converted_rankscore","SIFT_pred","Uniprot_acc_Polyphen2",
"Uniprot_id_Polyphen2","Uniprot_aapos_Polyphen2","Polyphen2_HDIV_score",
"Polyphen2_HDIV_rankscore","Polyphen2_HDIV_pred","Polyphen2_HVAR_score",
"Polyphen2_HVAR_rankscore","Polyphen2_HVAR_pred","LRT_score","LRT_converted_rankscore",
"LRT_pred","LRT_Omega","MutationTaster_score","MutationTaster_converted_rankscore",
"MutationTaster_pred","MutationTaster_model","MutationTaster_AAE",
"MutationAssessor_UniprotID","MutationAssessor_variant",
"MutationAssessor_score","MutationAssessor_score_rankscore",
"MutationAssessor_pred","FATHMM_score","FATHMM_converted_rankscore",
"FATHMM_pred","PROVEAN_score","PROVEAN_converted_rankscore","PROVEAN_pred",
"Transcript_id_VEST3","Transcript_var_VEST3","VEST3_score","VEST3_rankscore",
"MetaSVM_score","MetaSVM_rankscore","MetaSVM_pred","MetaLR_score","MetaLR_rankscore",
"MetaLR_pred","Reliability_index","M-CAP_score","M-CAP_rankscore","M-CAP_pred",
"REVEL_score","REVEL_rankscore","MutPred_score","MutPred_rankscore","MutPred_protID",
"MutPred_AAchange","MutPred_Top5features","CADD_raw","CADD_raw_rankscore",
"CADD_phred","DANN_score","DANN_rankscore","fathmm-MKL_coding_score",
"fathmm-MKL_coding_rankscore","fathmm-MKL_coding_pred",
"fathmm-MKL_coding_group","Eigen_coding_or_noncoding",
"Eigen-raw","Eigen-phred","Eigen-PC-raw","Eigen-PC-phred",
"Eigen-PC-raw_rankscore","GenoCanyon_score","GenoCanyon_score_rankscore",
"integrated_fitCons_score","integrated_fitCons_score_rankscore",
"integrated_confidence_value","GM12878_fitCons_score",
"GM12878_fitCons_score_rankscore","GM12878_confidence_value",
"H1-hESC_fitCons_score","H1-hESC_fitCons_score_rankscore",
"H1-hESC_confidence_value","HUVEC_fitCons_score",
"HUVEC_fitCons_score_rankscore","HUVEC_confidence_value",
"GERP++_NR","GERP++_RS","GERP++_RS_rankscore",
"phyloP100way_vertebrate","phyloP100way_vertebrate_rankscore",
"phyloP20way_mammalian","phyloP20way_mammalian_rankscore",
"phastCons100way_vertebrate","phastCons100way_vertebrate_rankscore",
"phastCons20way_mammalian","phastCons20way_mammalian_rankscore",
"SiPhy_29way_pi","SiPhy_29way_logOdds","SiPhy_29way_logOdds_rankscore",
"1000Gp3_AC","1000Gp3_AF","1000Gp3_AFR_AC","1000Gp3_AFR_AF",
"1000Gp3_EUR_AC","1000Gp3_EUR_AF","1000Gp3_AMR_AC","1000Gp3_AMR_AF",
"1000Gp3_EAS_AC","1000Gp3_EAS_AF","1000Gp3_SAS_AC","1000Gp3_SAS_AF",
"TWINSUK_AC","TWINSUK_AF","ALSPAC_AC","ALSPAC_AF","ESP6500_AA_AC",
"ESP6500_AA_AF","ESP6500_EA_AC","ESP6500_EA_AF","ExAC_AC","ExAC_AF",
"ExAC_Adj_AC","ExAC_Adj_AF","ExAC_AFR_AC","ExAC_AFR_AF","ExAC_AMR_AC",
"ExAC_AMR_AF","ExAC_EAS_AC","ExAC_EAS_AF","ExAC_FIN_AC","ExAC_FIN_AF",
"ExAC_NFE_AC","ExAC_NFE_AF","ExAC_SAS_AC","ExAC_SAS_AF","ExAC_nonTCGA_AC",
"ExAC_nonTCGA_AF","ExAC_nonTCGA_Adj_AC","ExAC_nonTCGA_Adj_AF",
"ExAC_nonTCGA_AFR_AC","ExAC_nonTCGA_AFR_AF","ExAC_nonTCGA_AMR_AC",
"ExAC_nonTCGA_AMR_AF","ExAC_nonTCGA_EAS_AC","ExAC_nonTCGA_EAS_AF",
"ExAC_nonTCGA_FIN_AC","ExAC_nonTCGA_FIN_AF","ExAC_nonTCGA_NFE_AC",
"ExAC_nonTCGA_NFE_AF","ExAC_nonTCGA_SAS_AC","ExAC_nonTCGA_SAS_AF",
"ExAC_nonpsych_AC","ExAC_nonpsych_AF","ExAC_nonpsych_Adj_AC",
"ExAC_nonpsych_Adj_AF","ExAC_nonpsych_AFR_AC","ExAC_nonpsych_AFR_AF",
"ExAC_nonpsych_AMR_AC","ExAC_nonpsych_AMR_AF","ExAC_nonpsych_EAS_AC",
"ExAC_nonpsych_EAS_AF","ExAC_nonpsych_FIN_AC","ExAC_nonpsych_FIN_AF",
"ExAC_nonpsych_NFE_AC","ExAC_nonpsych_NFE_AF","ExAC_nonpsych_SAS_AC",
"ExAC_nonpsych_SAS_AF","gnomAD_exomes_AC","gnomAD_exomes_AN",
"gnomAD_exomes_AF","gnomAD_exomes_AFR_AC","gnomAD_exomes_AFR_AN",
"gnomAD_exomes_AFR_AF","gnomAD_exomes_AMR_AC","gnomAD_exomes_AMR_AN",
"gnomAD_exomes_AMR_AF","gnomAD_exomes_ASJ_AC","gnomAD_exomes_ASJ_AN",
"gnomAD_exomes_ASJ_AF","gnomAD_exomes_EAS_AC","gnomAD_exomes_EAS_AN",
"gnomAD_exomes_EAS_AF","gnomAD_exomes_FIN_AC","gnomAD_exomes_FIN_AN",
"gnomAD_exomes_FIN_AF","gnomAD_exomes_NFE_AC","gnomAD_exomes_NFE_AN",
"gnomAD_exomes_NFE_AF","gnomAD_exomes_SAS_AC","gnomAD_exomes_SAS_AN",
"gnomAD_exomes_SAS_AF","gnomAD_exomes_OTH_AC","gnomAD_exomes_OTH_AN",
"gnomAD_exomes_OTH_AF","gnomAD_genomes_AC","gnomAD_genomes_AN","gnomAD_genomes_AF",
"gnomAD_genomes_AFR_AC","gnomAD_genomes_AFR_AN","gnomAD_genomes_AFR_AF",
"gnomAD_genomes_AMR_AC","gnomAD_genomes_AMR_AN","gnomAD_genomes_AMR_AF",
"gnomAD_genomes_ASJ_AC","gnomAD_genomes_ASJ_AN","gnomAD_genomes_ASJ_AF",
"gnomAD_genomes_EAS_AC","gnomAD_genomes_EAS_AN","gnomAD_genomes_EAS_AF",
"gnomAD_genomes_FIN_AC","gnomAD_genomes_FIN_AN","gnomAD_genomes_FIN_AF",
"gnomAD_genomes_NFE_AC","gnomAD_genomes_NFE_AN","gnomAD_genomes_NFE_AF",
"gnomAD_genomes_OTH_AC","gnomAD_genomes_OTH_AN","gnomAD_genomes_OTH_AF",
"clinvar_rs","clinvar_clnsig","clinvar_trait","clinvar_golden_stars",
"Interpro_domain","GTEx_V6p_gene","GTEx_V6p_tissue"]


def gene_to_interval(gene, release=75):
    data = pyensembl.EnsemblRelease(release)
    ret = None
    try:
        ret = data.genes_by_name(gene)[0]
    except:
        sys.stderr.write(" ".join(["Gene not found: ", gene, "\n"]))
    if ret is not None:
        return (str(ret.contig), str(ret.start), str(ret.end), str(ret.name))
    else:
        return ret

def genes_to_intervals(gene_list):
    gti = gene_to_interval
    ret = []
    for i in gene_list:
        r = gti(i)
        if r is not None:
            ret.append(r)
    return ret

def q_tabix(filename, chrom, start, end):
    """Call tabix and generate an array of strings for each line it returns."""
    query = '{}:{}-{}'.format(chrom, start, end)
    process = subprocess.Popen(['tabix', '-f', filename, query], stdout=subprocess.PIPE)
    for line in process.stdout:
        yield line.strip().split()

## Present a numerical score as well as a thresholded
## effect prediciton for each variant
def score_variant(var, annos):
    return 2.0

def match_var(tokens, annotokens):
    #print "call", tokens[0:5]
    #print "anno", annotokens[0:5]
    
    if tokens[0] == annotokens[0] and \
    tokens[1] == annotokens[1] and \
    tokens[3] == annotokens[2] and \
    tokens[4] == annotokens[3]:
        return True
    return False

def handle_line(line):

    
    ## Check if we have received a header line
    if line.startswith("#"):
        return line.strip()
    
    ## Strip / split our line
    tokens = line.strip().split("\t")

    ## vcf filter
    if not (tokens[6] == "PASS" or tokens[6] == "."):
        sys.stderr.write("Fails filter: " + tokens[0] + " " + tokens[1] + "\n") 
        return None

    ## Query tabix db
    tb_annos = None
    try:
        #tb_annos = tb_db.query(tokens[0], int(tokens[1]) - 1, int(tokens[1]) - 1 + max(len(tokens[4]), len(tokens[5])))
        tb_annos = q_tabix(default_tb_database, tokens[0], int(tokens[1]) - 1, int(tokens[1]) - 1 + max(len(tokens[4]), len(tokens[5])))
    except Exception:
        print "Error accessing tabix"
        pass

    modded = False
    ## ANNOTATE
    if tb_annos is not None:
        count = 0
        for i in tb_annos:
            #tb_toks = i.split("\t")
            tb_toks = i
            if match_var(tokens, tb_toks):
                sys.stderr.write("Variant found in DB:" + tokens[0] + " " + tokens[1] + " " + tokens[3] + "->" + tokens[4] + "\n")
                tinfo = []
                ## Make a bunch of infos and append to info string
                for j in xrange(0, len(header)):
                    if tb_toks[j] != ".":
                        tinfo.append("=".join([header[j], tb_toks[j]]))
                tokens[7] = tokens[7] + ";" + ";".join(tinfo)
                modded = True

    ## return annotated variant or None

    return "\t".join(tokens)

def bcsq(vfile, isPhased=False):
    cmd = "bcftools csq" + "-f " + fasta_ref
    if isPhased:
        cmd += " -l "
        cmd + " -g " + gtf_ref + " " + vcfile + " > tmp && mv " + vfile
    return vfile

def sort_vcf(vfi, faidx = None):
    oname = ".".join("".join(vfi.split("/")[-1]).split(".")[:-1]) + ".sorted.vcf"
    cmd = "bedtools sort -header "
    if faidx is not None:
        cmd += " -faidx " + faidx
    cmd += " -i " + vfi + "> " + oname
    sys.stderr.write("Sorting VCF file...\n")
    subprocess.call(cmd, shell = True)
    sys.stderr.write("Done.\n")
    return oname

def file_handler(fi):
    with open(fi, "r") as ifi:
        for line in ifi:
            yield line

def p_annotate_vars(vfi, runSerial = False):
    
    ## The serial appraoch chains generators
    ## I don't know if this is legal or advisable
    ## TODO subject to change.
    if runSerial:
        for i in file_handler(vfi):
            print handle_line(i)
        return []
    else:
        pool = multiprocessing.Pool(4)
        hl = handle_line
        lines = pool.map(hl, file_handler(vfi))
        return lines

if __name__ == "__main__":
    
    if args.tabix is not None:
        global default_tb_database
        default_tb_database = args.tabix.name
        tb_db = tabix.open(args.tabix.name)
    else:
        tb_db = tabix.open(default_tb_database)

    if args.gene_list is not None:
        pass
    else:
        pass

    if args.bedfile is not None:
        pass
    else:
        pass

    if args.infile is not None:
        vfi = sort_vcf(args.infile.name)
        lines = p_annotate_vars(vfi)
        for i in lines:
            if i is not None:
                print i
    else:
        sys.stderr.write("No input file. Please provide one (a bgzipped/tabixed VCF) with -i\n")
        exit(1)
