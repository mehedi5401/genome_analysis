#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH -J deseq2
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

module load R/4.5

mkdir -p /home/mdha5802/R/library

Rscript /home/mdha5802/genome_analysis/code/differential_expression/rnaseq/differential_expression_rnaseq.R