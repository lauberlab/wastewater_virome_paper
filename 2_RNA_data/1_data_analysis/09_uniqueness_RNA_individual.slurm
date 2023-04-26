#!/bin/bash
#SBATCH --output=RNA/logs/uniqueness_individual.out
#SBATCH --error=RNA/logs/uniqueness_individual.err
#SBATCH --time=11:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=1970
#SBATCH --job-name=MDCuniqRNA
#SBATCH --partition=romeo

# directories
wdir="RNA/uniqueness_individual"
mkdir ${wdir}
cd    ${wdir}

# extract true viral sequences
#qu="viral.fasta"
#for sd in ../assemblies/*; do
#	q0="${sd}/mmseqs/viral.fasta"
#	t0="${sd}/mmseqs/mmseqs_taxonomy.tsv"
#	grep 'Archaea;'   ${t0} | awk '{print $1;}' | sed 's/_[[:digit:]]\+$//' | sort | uniq  > tmp.ids
#	grep 'Bacteria;'  ${t0} | awk '{print $1;}' | sed 's/_[[:digit:]]\+$//' | sort | uniq >> tmp.ids
#	grep 'Eukaryota;' ${t0} | awk '{print $1;}' | sed 's/_[[:digit:]]\+$//' | sort | uniq >> tmp.ids
#	extractFromFasta.pl ${q0} tmp.ids -invert >> ${qu}
#	rm -rf tmp.ids
#done

# remove sequence redundancy via clustering
#mmseqs easy-linclust ${qu} viral_clustered tmp_mmseqs --threads 24 --min-seq-id 0.95 -c 0.65 --cluster-mode 2
#rm -rf tmp_mmseqs
#mv     viral_clustered_rep_seq.fasta viral_unique.fasta0
#rm -rf viral_clustered_*
#../../scripts/rename_contigs.pl viral_unique.fasta0 > viral_unique.fasta
#rm -rf viral_unique.fasta0
qu="viral_unique.fasta"

# build bowtie2 index
#bowtie2-build ${qu} viral_unique

# homology searches
pa=" --ultra-sensitive --threads 24 --masking 0 -k 1 -f 6 qseqid qlen sseqid pident evalue cigar "
diamond blastx -q ${qu} -d /beegfs/p_sra/db/RefSeq/viral_protein.dmnd                      ${pa} > viralprotein_diamond.tsv
#diamond blastx -q ${qu} -d /beegfs/p_sra/db/mmseqs2/viral_protein/viral_protein_nr.dmnd    ${pa} > viralprotein_diamond.tsv
diamond blastx -q ${qu} -d /beegfs/p_sra/db/Serratus/rdrp1.dmnd                            ${pa} > rdrp1_diamond.tsv
diamond blastx -q ${qu} -d /beegfs/p_sra/db/Serratus/ambi_rdrp_id90.dmnd                   ${pa} > ambi_diamond.tsv
diamond blastx -q ${qu} -d /beegfs/p_sra/db/Serratus/palmdb_uniques.dmnd                   ${pa} > palmdb_diamond.tsv
diamond blastx -q ${qu} -d /beegfs/p_sra/db/TaraOceans_Sullivan/RdRp_footprints_Tara.dmnd  ${pa} > tara_diamond.tsv
diamond blastx -q ${qu} -d /beegfs/p_sra/db/RdRp-scan/RdRp-scan_0.90.dmnd                  ${pa} > rdrpscan_diamond.tsv
diamond blastx -q ${qu} -d /beegfs/p_sra/db/NVPC_Koonin/RdRPs_ali822x.dmnd                 ${pa} > rvmt_diamond.tsv

