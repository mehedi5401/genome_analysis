#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 040:20:00
#SBATCH -J trim_rna_qc
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out


# Load modules
module load FastQC/0.12
module load MultiQC/1.28

# ── Paths ────────────────────────────────────────────────────────
DATADIR=/proj/uppmax2026-1-61/nobackup/mdha5802/trimmed/transcriptomics
BHI=$DATADIR/BHI
SERUM=$DATADIR/Serum
OUTDIR=/home/mdha5802/genome_analysis/output/qc/transcriptomics_data/trimmed_rna_seq

mkdir -p $OUTDIR/BHI
mkdir -p $OUTDIR/Serum
mkdir -p $OUTDIR/multiqc

# ── FastQC on BHI raw reads ──────────────────────────────────────
# echo "Running FastQC on BHI samples..."
# fastqc $BHI/*.fastq.gz --outdir $OUTDIR/BHI --threads 4
# echo "BHI done"

# ── FastQC on Serum raw reads ────────────────────────────────────
echo "Running FastQC on Serum samples..."
fastqc $SERUM/*.fastq.gz --outdir $OUTDIR/Serum --threads 4
echo "Serum done"

echo "All FastQC complete"
echo "Reports in: $OUTDIR"

echo "Running MultiQC..."
multiqc \
    $OUTDIR/BHI \
    $OUTDIR/Serum \
    --outdir $OUTDIR/multiqc \
    --filename rnaseq_multiqc_report
echo "MultiQC done"
