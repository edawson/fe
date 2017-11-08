fe
----------
Eric T Dawson  
October 2017

### Introduction
fe is a single python script that annotates VCF variants with predicted functional impacts.
These come from an aggregated database (dbNSFP) containing variant effect predictions from
multiple tools and optionally from the raw output of bcftools csq (a fast, phase-aware VEP clone from Petr Danecek).
Bedtools (Aaron Quinlan + others)  is used to sort and intersect VCF files.

### Installation
fe requires wooey (a python web app framework), pytabix, pyensembl, bedtools and bcftools.


### Usage

        python scan.py -i my_vcf.vcf -t tabixFI.gz


### Getting help
post on [github](https://github.com/edawson/fe) for bugs and help!
