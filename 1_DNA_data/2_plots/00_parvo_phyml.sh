#!/bin/sh


aln="parvo_NS1_trim"
fafile1="data/${aln}.fasA"
fafile2="data/${aln}_rn.fasA"
phyfile="data/${aln}.phy"

sed 's/>ref__/>/g' $fafile1 > $fafile2

fasta2relaxedPhylip.pl -f $fafile2 -o $phyfile
phyml -i $phyfile -d aa -m LG -f m -v 0 -a e -c 4 -o tlr
mv ${phyfile}_phyml_stats.txt ${phyfile}_phyml_SH_stats.txt
mv ${phyfile}_phyml_tree.txt  ${phyfile}_phyml_SH_tree.txt

