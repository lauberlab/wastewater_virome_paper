#!/bin/bash
#SBATCH --array=1-1000
#SBATCH --output=DNA/logs/%a_hmmsearch.out
#SBATCH --error=DNA/logs/%a_hmmsearch.err
#SBATCH --time=7:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=1900
#SBATCH --job-name=MDChmmsearch

# INPUT
D="DNA" # DNA or RNA

# get order name
sf=$(grep -w ${SLURM_ARRAY_TASK_ID} ${D}/hmmsearch/split_files.txt | awk '{print $2}')

# run HMMER
hmmsearch --cpu 6 --tblout ${sf}_NCLDVs_hittable.tsv          /home/lauber/lauber-2015-virusHunter/profiles/NCLDVs_lineagesMCP_n74.hmm  ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_AsfarPox_hittable.tsv        /home/lauber/lauber-2015-virusHunter/profiles/Asfar.Pox_MCP.hmm           ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_smallDNAviruses_hittable.tsv /home/lauber/lauber-2015-virusHunter/profiles/smallDNAviruses-profile.hmm ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_AdenoHerpes_hittable.tsv     /home/lauber/lauber-2015-virusHunter/profiles/AdenoHerpes-n11.hmm         ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_VirSorter2_hittable.tsv      /beegfs/p_sra/db/VirSorter2/hmm/viral/combined.hmm                        ${sf} > /dev/null

