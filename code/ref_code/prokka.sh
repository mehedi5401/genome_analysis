#!/bin/bash
#SBATCH --job-name=prokka
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=prokka_%j.log

module load prokka
# or: conda activate prokka

mkdir -p 1_Zhang_2017/annotation/prokka_output

prokka \
    1_Zhang_2017/assembly/unicycler_output/assembly.fasta \
    --outdir 1_Zhang_2017/annotation/prokka_output \
    --genus Enterococcus \
    --species faecium \
    --strain E745 \
    --prefix E745 \
    --cpus 8 \
    --kingdom Bacteria

echo "Annotation complete"