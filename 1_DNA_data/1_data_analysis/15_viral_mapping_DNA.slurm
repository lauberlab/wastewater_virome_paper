#!/bin/bash
#SBATCH --array=1-24
#SBATCH --output=DNA/logs/%a_viral_mapping.out
#SBATCH --error=DNA/logs/%a_viral_mapping.err
#SBATCH --time=7:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=1970
#SBATCH --job-name=mapMDCdna

# get order name
dataset=$(grep -w ${SLURM_ARRAY_TASK_ID} datasets_DNA.txt | awk '{print $2}')

# directories and files
ddir="DNA/fastq"
wdir="DNA/mapping"
 ref="DNA/uniqueness/viral_unique"

# run bowtie2
#  -L : length of seed alignment = first part of read with highest quality bases (bowtie2 default: 22)
#  -N : number of allow mismatches in seed alignment (bowtie2 default: 0)
bowtie2 --threads 12 --no-unal -L 20 -N 1 -x ${ref} -1 ${ddir}/${dataset}_R1_trim.fastq.gz -2 ${ddir}/${dataset}_R2_trim.fastq.gz -S ${wdir}/${dataset}.sam
samtools sort -@ 12 -o ${wdir}/${dataset}.bam ${wdir}/${dataset}.sam
samtools index         ${wdir}/${dataset}.bam
rm -rf                 ${wdir}/${dataset}.sam

# count reads
samtools idxstats ${wdir}/${dataset}.bam | cut -f 1,3 > ${wdir}/${dataset}.counts

