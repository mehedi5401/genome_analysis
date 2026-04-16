#!/bin/bash -l
#SBATCH -A uppmax2026-1-61

# Load modules
module load FastQC/0.12
module load MultiQC/1.28

# ── Paths ────────────────────────────────────────────────────────
OUTPUT=/home/mdha5802/genome_analysis/output/qc
DATADIR=/home/mdha5802/1_Zhang_2017
BHI=$DATADIR/transcriptomics_data/RNA-Seq_BH/raw
SERUM=$DATADIR/transcriptomics_data/RNA-Seq_Serum/raw
OUTDIR=$OUTPUT/transcriptomics_data

mkdir -p $OUTDIR/BHI
mkdir -p $OUTDIR/Serum
mkdir -p $OUTDIR/multiqc

# ── FastQC on BHI raw reads ──────────────────────────────────────
# echo "Running FastQC on BHI samples..."
# fastqc \
#     $BHI/ERR1797972_1.fastq.gz \
#     $BHI/ERR1797972_2.fastq.gz \
#     $BHI/ERR1797973_1.fastq.gz \
#     $BHI/ERR1797973_2.fastq.gz \
#     $BHI/ERR1797974_1.fastq.gz \
#     $BHI/ERR1797974_2.fastq.gz \
#     --outdir $OUTDIR/BHI \
#     --threads 8
# echo "BHI done"

# ── FastQC on Serum raw reads ────────────────────────────────────
# echo "Running FastQC on Serum samples..."
# fastqc \
#     $SERUM/ERR1797969_1.fastq.gz \
#     $SERUM/ERR1797969_2.fastq.gz \
#     $SERUM/ERR1797970_1.fastq.gz \
#     $SERUM/ERR1797970_2.fastq.gz \
#     $SERUM/ERR1797971_1.fastq.gz \
#     $SERUM/ERR1797971_2.fastq.gz \
#     --outdir $OUTDIR/Serum \
#     --threads 8
# echo "Serum done"

echo "All FastQC complete"
echo "Reports in: $OUTDIR"

echo "Running MultiQC..."
multiqc \
    $OUTDIR/BHI \
    $OUTDIR/Serum \
    --outdir $OUTDIR/multiqc \
    --filename rnaseq_multiqc_report
echo "MultiQC done"
