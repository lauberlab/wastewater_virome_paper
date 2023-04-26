#!/bin/bash
#SBATCH --output=DNA/logs/spades.out
#SBATCH --error=DNA/logs/spades.err
#SBATCH --time=95:59:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=54000
#SBATCH --job-name=MDCspades
#SBATCH --partition=julia

# run assembly script
03_joint_assembly_DNA.pl DNA 64 54 1

