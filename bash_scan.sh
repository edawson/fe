
bedtools=bedtools
tabix=tabix
faidx=~/Dropbox/refs/Homo_sapiens_assembly19.fasta.fai
annos=~/sandbox/fe/dbNSFP3.5a_variant.chr1.gz

vcf=$1

bedtools sort -header -faidx $faidx -i $vcf > $(basename $vcf .vcf).sorted.vcf && \
python vcf_filt.py $(basename $vcf .vcf).sorted.vcf > $(basename $vcf .vcf).sorted.filtered.vcf && \
for i in `grep -v "#" $(basename $vcf .vcf).sorted.filtered.vcf | cut -f 1,2 | gsed "s/\t/_/g"`
do 
    tabix $annos $(echo $i | cut -f 1 -d "_"):$(echo $i | cut -f 2 -d "_")-$(echo $i | cut -f 2 -d "_") >> annos.txt
done

