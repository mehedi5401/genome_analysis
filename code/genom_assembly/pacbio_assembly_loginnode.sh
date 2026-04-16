#!/bin/bash -l
#SBATCH -A uppmax2026-1-61


# Load modules
module load canu/2.3
module load SAMtools/1.22

# Define paths
WORKDIR=/home/mdha5802/1_Zhang_2017
PACBIO=$WORKDIR/genomics_data/PacBio
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/assembly/pacbio

# Create output directory
mkdir -p $OUTDIR

# Step 1 — Combine all PacBio subread files
echo "Combining PacBio files..."
cat $PACBIO/m131023*.fastq.gz \
    $PACBIO/m131024*.fastq.gz \
    > $OUTDIR/pacbio_combined.fastq.gz
echo "Combined size: $(du -sh $OUTDIR/pacbio_combined.fastq.gz | cut -f1)"

# Step 2 — Run Canu assembly
echo "Starting Canu assembly..."
canu \
    -p E745 \
    -d $OUTDIR \
    genomeSize=2.8m \
    minReadLength=1000 \
    minOverlapLength=500 \
    useGrid=false\
    -pacbio-raw $OUTDIR/pacbio_combined.fastq.gz


echo "Canu assembly complete"
echo "Output files:"
ls -lh $OUTDIR