#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH -J flye_assembly
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools
module load Flye

# Define paths
WORKDIR=/home/mdha5802/genome_analysis/data/PacBio
PACBIO=$WORKDIR
OUTDIR=$WORKDIR

# Create output directory
# mkdir -p $OUTDIR

# Step 1 — Combine PacBio subread files into one
# echo "Combining PacBio files..."
# cat $PACBIO/m131023*.fastq.gz $PACBIO/m131024*.fastq.gz > $PACBIO/pacbio_combined.fastq.gz
# echo "Done combining. Size: $(du -sh $PACBIO/pacbio_combined.fastq.gz | cut -f1)"

# Step 2 — Run Flye hybrid assembly
echo "Starting Flye assembly..."
flye \
    --pacbio-raw $PACBIO/pacbio_combined.fastq.gz \
    --out-dir $OUTDIR \
    --genome-size 2.8m \
    --threads 16

echo "Flye assembly complete"
echo "Output files:"
ls -lh $OUTDIR