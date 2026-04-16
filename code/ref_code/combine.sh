#!/bin/bash
#SBATCH --job-name=combine_pacbio
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=combine_pacbio_%j.log

# Combine all PacBio subreads into one file
cat 1_Zhang_2017/genomics_data/PacBio/*.fastq.gz > \
    1_Zhang_2017/genomics_data/PacBio/pacbio_combined.fastq.gz

echo "Done combining PacBio files"