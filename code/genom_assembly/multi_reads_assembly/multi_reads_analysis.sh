#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH -J spades_assembly
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

# Load modules
module load SPAdes/4.2

# ── Paths ────────────────────────────────────────────────────────
DATADIR=/home/mdha5802/1_Zhang_2017
ILLUMINA=$DATADIR/genomics_data/Illumina
NANO=$DATADIR/genomics_data/Nanopore
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/assembly/illumina_nano

mkdir -p $OUTDIR

echo "Starting SPAdes assembly..."

spades.py \
    --pe1-1 $ILLUMINA/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
    --pe1-2 $ILLUMINA/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
    --nanopore $NANO/E745_all.fasta.gz \
    -o $OUTDIR \
    --threads 4 \
    --memory 32

echo "SPAdes assembly complete"
echo "Output files:"
ls -lh $OUTDIR