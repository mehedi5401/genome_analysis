#!/bin/bash
#SBATCH --job-name=unicycler_assembly
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=unicycler_%j.log

# Load unicycler - exact command depends on your HPC
# Common options:
module load unicycler
# or: conda activate unicycler

# Create output directory
mkdir -p 1_Zhang_2017/assembly/unicycler_output

# Run hybrid assembly
unicycler \
    -1 1_Zhang_2017/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
    -2 1_Zhang_2017/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
    -l 1_Zhang_2017/genomics_data/Nanopore/E745_all.fasta.gz \
    --unpaired 1_Zhang_2017/genomics_data/PacBio/pacbio_combined.fastq.gz \
    -o 1_Zhang_2017/assembly/unicycler_output \
    -t 16

echo "Assembly complete"