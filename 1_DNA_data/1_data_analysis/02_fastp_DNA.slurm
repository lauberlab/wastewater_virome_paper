#!/bin/bash
#SBATCH --array=1-24
#SBATCH --output=DNA/logs/%a_fastp.out
#SBATCH --error=DNA/logs/%a_fastp.err
#SBATCH --time=0:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=1800
#SBATCH --job-name=MDCfastp

# INPUT
D="DNA" # DNA or RNA

# get order name
dataset=$(grep -w ${SLURM_ARRAY_TASK_ID} datasets_${D}.txt | awk '{print $2}')

# files and directories
name="ww_${D}_all"
data="/beegfs/p_sra/data/MDC_wastewater/${name}.tar"
wdir=`pwd`

# extract data set
mkdir tmp/${dataset}
tar -xvf $data -C tmp/${dataset}/ ${name}/${dataset}_R1_001.fastq.gz
tar -xvf $data -C tmp/${dataset}/ ${name}/${dataset}_R2_001.fastq.gz

# run fastp
cd   tmp/${dataset}/${name}
fastp -i ${dataset}_R1_001.fastq.gz  -I ${dataset}_R2_001.fastq.gz  \
      -o ${dataset}_R1_trim.fastq.gz -O ${dataset}_R2_trim.fastq.gz \
      -q 20 -w 12 --dedup
cd   $wdir
mv   tmp/${dataset}/${name}/*_trim.fastq.gz ${D}/fastq/
rm   -rf tmp/${dataset}

