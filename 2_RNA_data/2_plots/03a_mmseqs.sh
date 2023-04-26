#!/bin/sh

parameters1="-c 0.8 -e 1e-5 --min-seq-id 0.95"
mmseqs easy-cluster data/lenar_known_and_new_palmscan.fasta   clustered_0 tmp_mmseqs $parameters1 --threads 4
rm -rf tmp_mmseqs

parameters2="-c 0.8 -e 1e-5 --min-seq-id 0.5"
mmseqs easy-cluster clustered_0_rep_seq.fasta                 clustered_1 tmp_mmseqs $parameters2 --threads 4
rm -rf tmp_mmseqs

