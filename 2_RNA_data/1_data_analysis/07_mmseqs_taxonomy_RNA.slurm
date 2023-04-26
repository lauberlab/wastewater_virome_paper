#!/bin/bash
#SBATCH --output=RNA/logs/mmseqs_taxonomy_unique.out
#SBATCH --error=RNA/logs/mmseqs_taxonomy_unique.err
#SBATCH --time=71:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10000
#SBATCH --job-name=MDCtaxRNA
#SBATCH --partition=haswell256

# directories
wdir="RNA/mmseqs_unique"

# files
infile="RNA/uniqueness_individual/viral_unique.fasta"

# extract protein sequences
getorf -sequence $infile -outseq ${wdir}/viral_unique_aa.fasta -minsize 300

# mmseqs taxonomy
cd ${wdir}
mmseqsdb="/beegfs/p_sra/db/mmseqs2/nr/nr.fa.db"
mmseqs createdb       viral_unique_aa.fasta viral_unique_aa.db --dbtype 1
mmseqs taxonomy       viral_unique_aa.db $mmseqsdb mmseqs_taxonomy_unique tmp_taxonomy --threads 56 --lca-ranks superkingdom,kingdom,phylum,subphylum,class,order,suborder,family,subfamily,genus,subgenus,species
mmseqs createtsv      viral_unique_aa.db mmseqs_taxonomy_unique mmseqs_taxonomy_unique.tsv
mmseqs taxonomyreport $mmseqsdb   mmseqs_taxonomy_unique mmseqs_taxonomy_unique_for_pavian.txt
rm -rf tmp_taxonomy

