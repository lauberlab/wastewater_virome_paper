#!/bin/bash
#SBATCH --output=DNA/logs/mmseqs_taxonomy.out
#SBATCH --error=DNA/logs/mmseqs_taxonomy.err
#SBATCH --time=95:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=54000
#SBATCH --job-name=MDCtaxDNA
#SBATCH --partition=julia

# directories
wdir="DNA/mmseqs"
mkdir ${wdir}
cd    ${wdir}

# combine results of homology searches
head -n 3 ../hmmsearch/split/getorf.fasta.1_NCLDVs_hittable.tsv          > ../hmmsearch/NCLDVs_hittable.tsv
head -n 3 ../hmmsearch/split/getorf.fasta.1_AsfarPox_hittable.tsv        > ../hmmsearch/AsfarPox_hittable.tsv
head -n 3 ../hmmsearch/split/getorf.fasta.1_smallDNAviruses_hittable.tsv > ../hmmsearch/smallDNAviruses_hittable.tsv
head -n 3 ../hmmsearch/split/getorf.fasta.1_AdenoHerpes_hittable.tsv     > ../hmmsearch/AdenoHerpes_hittable.tsv
head -n 3 ../hmmsearch/split/getorf.fasta.1_VirSorter2_hittable.tsv      > ../hmmsearch/VirSorter2_hittable.tsv

cat ../hmmsearch/split/*_NCLDVs_hittable.tsv          | grep -v '#' >> ../hmmsearch/NCLDVs_hittable.tsv
cat ../hmmsearch/split/*_AsfarPox_hittable.tsv        | grep -v '#' >> ../hmmsearch/AsfarPox_hittable.tsv
cat ../hmmsearch/split/*_smallDNAviruses_hittable.tsv | grep -v '#' >> ../hmmsearch/smallDNAviruses_hittable.tsv
cat ../hmmsearch/split/*_AdenoHerpes_hittable.tsv     | grep -v '#' >> ../hmmsearch/AdenoHerpes_hittable.tsv
cat ../hmmsearch/split/*_VirSorter2_hittable.tsv      | grep -v '#' >> ../hmmsearch/VirSorter2_hittable.tsv

# combine hit lists and extract
cat ../hmmsearch/*_hittable.tsv | grep -v '#' | awk '{print $1;}' | sed 's/_[[:digit:]]\+$//' | sort | uniq > viral.ids
cat ../hmmsearch/*_hittable.tsv | grep -v '#' | awk '{print $1;}' |                             sort | uniq > viral_aa.ids
extractFromFasta.pl ../spades/scaffolds.fasta viral.ids    > viral.fasta
extractFromFasta.pl ../hmmsearch/getorf.fasta viral_aa.ids > viral_aa.fasta
rm viral.ids viral_aa.ids

# mmseqs taxonomy
mmseqsdb="/beegfs/p_sra/db/mmseqs2/nr/nr.fa.db"
mmseqs createdb       viral_aa.fasta viral_aa.db --dbtype 1
mmseqs taxonomy       viral_aa.db $mmseqsdb mmseqs_taxonomy tmp_taxonomy --threads 64 --lca-ranks superkingdom,kingdom,phylum,subphylum,class,order,suborder,family,subfamily,genus,subgenus,species
mmseqs createtsv      viral_aa.db mmseqs_taxonomy mmseqs_taxonomy.tsv
mmseqs taxonomyreport $mmseqsdb   mmseqs_taxonomy mmseqs_taxonomy_for_pavian.txt
rm -rf tmp_taxonomy

