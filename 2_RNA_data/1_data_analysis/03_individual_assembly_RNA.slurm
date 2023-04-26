#!/bin/bash
#SBATCH --array=1-116
#SBATCH --output=RNA/logs/%a_assembly.out
#SBATCH --error=RNA/logs/%a_assembly.err
#SBATCH --time=11:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=54000
#SBATCH --job-name=assMDCrna
#SBATCH --partition=julia

# get order name
dataset=$(grep -w ${SLURM_ARRAY_TASK_ID} datasets_RNA.txt | awk '{print $2}')

# directories
ddir="/beegfs/p_sra/lauber-2015-virusGatherer/14_MDC_wastewater_untargeted/RNA/fastq"
wdir="RNA/assemblies/${dataset}"
mkdir ${wdir}
cd    ${wdir}

# run spades
spades --checkpoints last -1 ${ddir}/${dataset}_R1_trim.fastq.gz -2 ${ddir}/${dataset}_R2_trim.fastq.gz -t 6 -m 324 -o spades

