#!/bin/bash
#SBATCH --job-name=fastqc_illumina
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=fastqc_illumina_%j.log

module load fastqc

mkdir -p 1_Zhang_2017/QC/fastqc/genomics/illumina_raw

fastqc \
    1_Zhang_2017/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
    1_Zhang_2017/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
    --outdir 1_Zhang_2017/QC/fastqc/genomics/illumina_raw \
    --threads 4

echo "FastQC Illumina genomic done"