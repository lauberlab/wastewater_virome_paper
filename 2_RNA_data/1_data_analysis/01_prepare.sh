#!/bin/sh

# input
infile1="/beegfs/p_sra/data/MDC_wastewater/ww_DNA_all.tar"
infile2="/beegfs/p_sra/data/MDC_wastewater/ww_RNA_all.tar"

# directories
mkdir DNA
mkdir DNA/logs
mkdir DNA/assemblies
mkdir RNA
mkdir RNA/logs
mkdir RNA/assemblies

# data sets
tar -tvf $infile1 | awk '{print $6;}' | sed 's/ww_DNA_all\///g' | sed 's/_R1_001.fastq.gz//g' | sed 's/_R2_001.fastq.gz//g' \
	| sort | uniq | grep -v -e '^$' | awk '{print NR,$1;}' > datasets_DNA.txt
tar -tvf $infile2 | awk '{print $6;}' | sed 's/ww_RNA_all\///g' | sed 's/_R1_001.fastq.gz//g' | sed 's/_R2_001.fastq.gz//g' \
	| sort | uniq | grep -v -e '^$' | awk '{print NR,$1;}' > datasets_RNA.txt

