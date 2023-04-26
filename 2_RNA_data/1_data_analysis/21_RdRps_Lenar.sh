#!/bin/sh

# directories and files
idir="RNA/uniqueness_individual"
odir="../RdRp"
onam="rdrp_contigs"
rnam="lenar"

# only look at contigs
#  - with RdRp hits
#  - of length >1000 nt
cd $idir
cat ambi_diamond.tsv palmdb_diamond.tsv rdrp1_diamond.tsv rdrpscan_diamond.tsv tara_diamond.tsv \
	| awk '{ if( $2>999 ){ print $1;} }' | sort | uniq > ${odir}/${onam}.ids
helper/extractFromFasta.pl viral_unique.fasta ${odir}/rdrp_contigs.ids > ${odir}/${onam}.fasta


# group into viral phyla
threads=6
profiles="helper/RNAvirome.branch12345.representatives.5newPhylaOceanVirome.hmm"
cd $odir
helper/domainFromContigs.pl ${onam}.fasta ${profiles} ${onam}_sorted 1000

# run palmscan on sorted proteins
fasfile1="${onam}_sorted/2_grouped/Lenarviricota_RdRp/protein.fas"
fasfile2="${onam}_sorted/2_grouped/Lenarviricota_RdRp/nucleotide.fas"
palmscan -search_pp ${fasfile1} -rdrp -ppout ${onam}_${rnam}_palmscan.fa -report ${onam}_${rnam}_palmscan.txt -fevout ${onam}_${rnam}_palmscan.fev -threads ${threads}
palmscan -search_pp ${fasfile2} -rdrp -ppout ${onam}_${rnam}_palmscan.fa -report ${onam}_${rnam}_palmscan.txt -fevout ${onam}_${rnam}_palmscan.fev -threads ${threads}

# keep canonical and permuted RdRp found by palscan
grep 'order=ABC' ${onam}_${rnam}_palmscan.fev | awk '{print $2;}' | awk 'BEGIN{FS="=";}{print $2;}' > ${rnam}_canonical_ABC.ids
helper/extractFromFasta.pl ${fasfile1} ${rnam}_canonical_ABC.ids > ${rnam}_canonical_ABC_protein.fasta
helper/extractFromFasta.pl ${fasfile2} ${rnam}_canonical_ABC.ids > ${rnam}_canonical_ABC_nucleotide.fasta

#grep 'order=CAB' ${onam}_${rnam}_palmscan.fev | awk '{print $2;}' | awk 'BEGIN{FS="=";}{print $2;}' > ${rnam}_permuted_CAB.ids
helper/extractFromFasta.pl ${fasfile1} ${rnam}_permuted_CAB.ids > ${rnam}_permuted_CAB_protein.fasta
helper/extractFromFasta.pl ${fasfile2} ${rnam}_permuted_CAB.ids > ${rnam}_permuted_CAB_nucleotide.fasta

# cluster
mmseqs easy-cluster ${rnam}_canonical_ABC_protein.fasta ${rnam}_canonical_ABC_protein tmp_mmseqs --min-seq-id 0.9 -c 0.65 --cluster-mode 2 --threads $threads
rm -rf tmp_mmseqs

