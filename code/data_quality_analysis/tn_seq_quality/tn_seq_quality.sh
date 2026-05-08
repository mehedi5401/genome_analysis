#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH -J tn_seq_qc
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out


# Load modules
module load FastQC/0.12
module load MultiQC/1.28

# ── Paths ────────────────────────────────────────────────────────
DATADIR=/home/mdha5802/1_Zhang_2017/transcriptomics_data
BHI=$DATADIR/Tn-Seq_BHI
HSERUM=$DATADIR/Tn-Seq_HSerum
SERUM=$DATADIR/Tn-Seq_Serum
OUTDIR=/home/mdha5802/genome_analysis/output/qc/transcriptomics_data/raw_tn_seq

mkdir -p $OUTDIR/BHI
mkdir -p $OUTDIR/HSerum
mkdir -p $OUTDIR/Serum
mkdir -p $OUTDIR/multiqc

# ── FastQC on BHI trimmed reads ──────────────────────────────────
echo "Running FastQC on BHI samples..."
fastqc $BHI/trim_*.fastq.gz --outdir $OUTDIR/BHI --threads 4
echo "BHI done"

# ── FastQC on HSerum trimmed reads ───────────────────────────────
echo "Running FastQC on HSerum samples..."
fastqc $HSERUM/trim_*.fastq.gz --outdir $OUTDIR/HSerum --threads 4
echo "HSerum done"

# ── FastQC on Serum trimmed reads ────────────────────────────────
echo "Running FastQC on Serum samples..."
fastqc $SERUM/trim_*.fastq.gz --outdir $OUTDIR/Serum --threads 4
echo "Serum done"

echo "All FastQC complete"
echo "Reports in: $OUTDIR"

# ── MultiQC across all three conditions ──────────────────────────
echo "Running MultiQC..."
multiqc \
    $OUTDIR/BHI \
    $OUTDIR/HSerum \
    $OUTDIR/Serum \
    --outdir $OUTDIR/multiqc \
    --filename tnseq_multiqc_report
echo "MultiQC done"