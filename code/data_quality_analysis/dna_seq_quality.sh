#!/bin/bash -l
#SBATCH -A uppmax2026-1-61

# Load modules
module load FastQC/0.12

# ── Paths ────────────────────────────────────────────────────────
DATADIR=/home/mdha5802/1_Zhang_2017/genomics_data/Illumina
OUTDIR=/home/mdha5802/genome_analysis/output/qc/genomics/raw_illumina

mkdir -p $OUTDIR

── FastQC on BHI raw reads ──────────────────────────────────────
echo "Running FastQC..."
fastqc \
    $DATADIR/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
    $DATADIR/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
    --outdir $OUTDIR \
    --threads 2
echo "Done"

echo "All FastQC complete"
echo "Reports in: $OUTDIR"
