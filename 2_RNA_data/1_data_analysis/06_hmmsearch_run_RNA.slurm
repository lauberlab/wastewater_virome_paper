#!/bin/bash
#SBATCH --array=1-4998
#SBATCH --output=RNA/logs/%a_hmmsearch.out
#SBATCH --error=RNA/logs/%a_hmmsearch.err
#SBATCH --time=7:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=1900
#SBATCH --job-name=MDChmmsearch

# INPUT
D="RNA" # DNA or RNA

# get order name
sf=$(grep -w ${SLURM_ARRAY_TASK_ID} ${D}/hmmsearch/split_files.txt | awk '{print $2}')

# run HMMER
hmmsearch --cpu 6 --tblout ${sf}_NCLDVs_hittable.tsv          /home/lauber/lauber-2015-virusHunter/profiles/NCLDVs_lineagesMCP_n74.hmm  ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_AsfarPox_hittable.tsv        /home/lauber/lauber-2015-virusHunter/profiles/Asfar.Pox_MCP.hmm           ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_smallDNAviruses_hittable.tsv /home/lauber/lauber-2015-virusHunter/profiles/smallDNAviruses-profile.hmm ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_AdenoHerpes_hittable.tsv     /home/lauber/lauber-2015-virusHunter/profiles/AdenoHerpes-n11.hmm         ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_VirSorter2_hittable.tsv      /beegfs/p_sra/db/VirSorter2/hmm/viral/combined.hmm                        ${sf} > /dev/null

hmmsearch --cpu 6 --tblout ${sf}_RNAviruses_hittable.tsv      /home/lauber/lauber-2015-virusHunter/profiles/RNAviruses-profile.hmm      ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_TaraOceans_hittable.tsv      /beegfs/p_sra/db/TaraOceans_Sullivan/RNAvirome_RdRp_5newPhyla.hmm         ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_NVPC_hittable.tsv            /beegfs/p_sra/db/NVPC_Koonin/NVPC.hmm                                     ${sf} > /dev/null
hmmsearch --cpu 6 --tblout ${sf}_RdRpScan_hittable.tsv        /beegfs/p_sra/db/RdRp-scan/RdRp-scan_0.90                                 ${sf} > /dev/null

