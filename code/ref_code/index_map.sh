#!/bin/bash
#SBATCH --job-name=bowtie2_index
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=bowtie2_index_%j.log

module load bowtie2

mkdir -p 1_Zhang_2017/rnaseq/genome_index

bowtie2-build \
    1_Zhang_2017/annotation/prokka_output/E745.fasta \
    1_Zhang_2017/rnaseq/genome_index/E745

echo "Index built"