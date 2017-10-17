import argparse
import subprocess
import os
import sys
from cStringIO import StringIO
import multiprocessing as mp
import tabix
from genelist_to_bed.genelist_to_bed import genes_to_intervals

dbSNP_file = ""
annotation_file = ""
ref_fasta = ""
faidx = "/Users/dawsonet/Dropbox/refs/Homo_sapiens_assembly19.fasta.fai"
ref_gtf = ""
bedtools = "bedtools"
bcftools = ""
bgzip = ""
#tabix = ""
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


#def parse_args():
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest = "infile", required = True, type = file)
parser.add_argument("-m", "--maf", dest = "maf", action = "store_true", help = "Input file is MAF format")
#parser.add_argument("-d", "--dbsnp", dest="dbsnp", action = "store_true", help = "Filter with dbSNP")
parser.add_argument("-t", "--tabix", dest="tabx", type = str, help = "A tabix-compatible file for annotation")
#parser.add_argument("-g", "--genome_build", dest = "genome", type = str, help = "The genome build to use")
parser.add_argument("-l", "--gene-list", dest="gene_list", type = str, help = "A file containing a specific list of genes to search")
parser.add_argument("-b", "--bed", dest="bed", action = "store_true", help = "Input is bed format")
parser.add_argument("-L", "--list-bed", dest = "list_bed", type = str, help = "A gene list bed file")
#parser.add_argument("-s", "--scan-bed", dest = "bedfilter", type = str, help = "A bedtools compatible file to filter variants against")
#parser.add_argument("-e", "--ensembl", dest = "ensembl", type = str, help = "The ensembl build to use")
parser.add_argument("-f", "--faidx", dest = "faidx", type = str, help = "FAIDX for vcf sorting")
#parser.add_argument("-r", "--reference", dest = "reference", type = str, help = "Reference genome")
#parser.add_argument("-p", "--phased", dest = "phased", action = "store_true", default = False, help = "VCF input is phased")
#    return parser.parse_args()

## Calculate an RMS of the sum of del annotations
def score_var(annos):
    score = 0.0
    for i in annos:
        score += float(i)
    return score

def list_to_bed(l):
    return "\n".join(["\t".join(i) for i in l])

## Take in an annotated VCF and append a composite score to each variant
def calculate_variant_scores(annotated_vcf):
    return

def escaped_bed_string(l):
    return "\n".join(["\t".join(i) for i in l])

def bedstring_select(query, bedstring):
    redirect_string = "<(echo -e " + "\"" + bedstring + "\")"
    #print redirect_string
    return bed_select(query, redirect_string)

def bed_filter(query, sieve):
    cmd = "bedtools intersect -header -v -a " + query + " -b " + sieve
    cap = subprocess.check_output(cmd)
    return cap

def tabix_annotate_variant_file(vfi, tfi):
    header = []
    annotated = []
    annos = []
    tb = tabix.open(tfi)
    with open(vfi, "r") as ifi:
        for line in ifi:
            if line.startswith("#"):
                header.append(line.strip())
                continue
            else:
                r = tabix_annotate_var(line, tb)
                annotated.append(r[0])
                if r[1] is not None:
                    annos.append(r[1])
    
    return header, annotated, annos


def annotate_variant_with_tabix(var, anno_tabix, cols=[0]):
    anno_d = {}
    for i in anno_tabix:
        if i[2] == var[3] and i[3] == var[4]:
            for c in cols:
                anno_d[header[c]] = i[c]
            if i[6] != ".":
                anno_d["dbSNP"] = i[6]
            if i[11] != ".":
                anno_d["Gene"] = i[11]
    n_info_l = []
    for i in anno_d:
        if anno_d[i] != ".":
            n_info_l.append(i)
    n_info_str = ";".join(["=".join([x, anno_d[x]]) for x in n_info_l])
    var[7] = var[7] + ";" + n_info_str
    print "\t".join(var)
    return var


def tabix_annotate_var(var, tb):
    i = var.split("\t")
    ret = None
    try:
        r = tb.query(i[0], int(i[1]), int(i[1]) + 1)
        ret = r
    except:
        pass
    if ret is not None:
        var = annotate_variant_with_tabix(i, ret, range(24, 32))
    return var, ret

def bed_select(query, sieve):
    oname = ".".join(query.split(".")[:-1]) + ".selected.vcf"
    cmd = "bedtools intersect -a " + query + " -b " + sieve + " > " + oname
    subprocess.call(cmd)
    return oname

def bcsq(vfile, isPhased=False):
    cmd = "bcftools csq" + "-f " + fasta_ref
    if isPhased:
        cmd += " -l "
    cmd + " -g " + gtf_ref + " " + vcfile + " > tmp && mv " + vfile
    return vfile

def sort_vcf(vfi, faidx):
    oname = ".".join("".join(vfi.split("/")[-1]).split(".")[:-1]) + ".sorted.vcf"
    cmd = "bedtools sort -header "
    if args.faidx is not None:
        cmd += " -faidx " + faidx
    cmd += " -i " + vfi + " > " + oname
    print cmd
    subprocess.call(cmd, shell = True)
    return oname

def pfunc(l):
    return l if (l.startswith("#") or line.split("\t")[6] in ["PASS", "."]) else None

def p_vcf_filter(vcf):
    tfi = open(vcf, "r")
    lines = tfi.readlines()
    lines = [pfunc(i) for i in lines]

def vcf_filter(vfi):
    ret = []
    with open(vfi, "r") as ifi:
        for line in ifi:
            line = line.strip()
            tokens = line.split("\t")
            if line.startswith("#") or tokens[6] == "PASS" or tokens[6] == ".":
                ret.append( line )
    ofname = ".".join(vfi.split(".")[:-1]) + ".filtered.vcf"
    with open(ofname, "w") as ofi:
        for r in ret:
            ofi.write(r + "\n")
    return ofname


def p_helper_func(vcf):
    for line in vcf:
        if not line.startswith("#"):
            yield line
## Return the header and an iterator (ok, a generator) to the main bit
def parse_vcf(vfi):
    header = []
    with open(vfi, "r") as ifi:
        for line in ifi:
            if line.startswith("#"):
                header.append(line)
    return header, p_helper_fun(vfi)

if __name__ == "__main__":
    args = parser.parse_args()

    infile = args.infile
   
    ## Sort whatever input we get and make it bedtools compatible
    ## MAF -> bed
    ## FASTQ will get mapped and called with BWA / FreeBayes
    variant_file = None
    if args.maf:
        print "Not implemented: MAF parsing"
        exit(9)
    else:
        variant_file = sort_vcf(args.infile, args.faidx)

    gene_bed = []
    gl = []
    ## Process any genes / beds we got, then filter our variants down to that set
    if args.gene_list is not None:
        with open(args.gene_list, "r") as ifi:
            for line in ifi:
                g = line.strip().split()[0]
                gl.append(g)
        gene_bed = genes_to_intervals(gl)
    if args.list_bed is not None:
        with open(args.list_bed, "r") as ifi:
            for line in ifi:
                line = line.strip()
                gene_bed.append(line.split("\t"))
                
    ## Filter vcf by gene_bed if it exists
    #if len(gene_bed) > 0:
    #    print list_to_bed(gene_bed)
    if len(gene_bed) > 1:
        variant_file = bedstring_select(variants, gene_bed)
    ## Filter VCF by VCF filter lines
    variant_file = vcf_filter(variant_file)

    #variant_file = pvcf_filter(variant_file)
    ## Run bcftools csq to predict mutation impacts at the amino acid level
    #variants = bcsq(variants)
    ## Pull annotations from annotation file and label VCF
    if args.tabx is not None:
        header, variants, annos = tabix_annotate_variant_file(variant_file, args.tabx)
        for x in annos:
            for y in x:
                print "\t".join(y)
    ## Score all variants based on annotations

    ## Produce a useful report
        




