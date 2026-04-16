#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=quast_%j.log

module load quast
# or: conda activate quast

mkdir -p 1_Zhang_2017/assembly/quast_output

quast.py \
    1_Zhang_2017/assembly/unicycler_output/assembly.fasta \
    -o 1_Zhang_2017/assembly/quast_output \
    --threads 4

echo "QUAST complete"