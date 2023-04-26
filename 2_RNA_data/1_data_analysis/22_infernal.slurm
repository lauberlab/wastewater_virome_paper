#!/bin/bash
#SBATCH --output=RNA/logs/infernal.out
#SBATCH --error=RNA/logs/infernal.err
#SBATCH --time=163:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=50000
#SBATCH --job-name=MDCinfernal
#SBATCH --partition=julia

# input
queryF="/beegfs/p_sra/db/Rfam/Rfam.cm"
ribozF="/beegfs/p_sra/db/Rfam/Rfam_ribozymes.txt"
infile1="RNA/RdRp/rdrp_contigs.fasta"
infile2="RNA/assemblies_combined/scaffolds.fasta"
outnam1="rdrp_contigs"
outnam2="all_contigs"
outdir="RNA/infernal"
threads=12

if [ ! -d $outdir ]; then
	mkdir $outdir
fi

# run infernal on RdRp contigs
cmsearch -E 0.001 --cpu ${threads} --tblout ${outdir}/${outnam1}_cmsearch.txt $queryF $infile1
grep -v '^#' ${outdir}/${outnam1}_cmsearch.txt | awk 'BEGIN{OFS="\t";}{print $1,$3,$4,$6,$7,$8,$9,$10,$15,$16;}' >  ${outdir}/${outnam1}_cmsearch.tsv
grep -w -f ${ribozF} ${outdir}/${outnam1}_cmsearch.tsv > ${outdir}/${outnam1}_cmsearch_ribozymes.tsv

# run infernal on all contigs
cmsearch -E 0.001 --cpu ${threads} --tblout ${outdir}/${outnam2}_cmsearch.txt $queryF $infile2 > /dev/null
grep -v '^#' ${outdir}/${outnam2}_cmsearch.txt | awk 'BEGIN{OFS="\t";}{print $1,$3,$4,$6,$7,$8,$9,$10,$15,$16;}' >  ${outdir}/${outnam2}_cmsearch.tsv
grep -w -f ${ribozF} ${outdir}/${outnam2}_cmsearch.tsv > ${outdir}/${outnam2}_cmsearch_ribozymes.tsv

