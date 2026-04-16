#!/bin/bash
#SBATCH --job-name=deseq2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=deseq2_%j.log

module load R

Rscript 1_Zhang_2017/rnaseq/deseq2_analysis.R