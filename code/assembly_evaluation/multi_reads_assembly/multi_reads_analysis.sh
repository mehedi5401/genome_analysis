#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 06:00:00
#SBATCH -J spades_assembly
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

# Load modules
module load SPAdes/4.2

# ── Paths ────────────────────────────────────────────────────────
DATADIR=/home/mdha5802/1_Zhang_2017
ILLUMINA=$DATADIR/genomics_data/Illumina
PACBIO=$DATADIR/genomics_data/PacBio
NANO=$DATADIR/genomics_data/Nanopore
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/assembly/multi_reads

mkdir -p $OUTDIR

echo "Starting SPAdes assembly..."

spades.py \
    --pe1-1 $ILLUMINA/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
    --pe1-2 $ILLUMINA/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
    --nanopore $NANO/E745_all.fasta.gz \
    --outdir $OUTDIR \
    --threads 4 \
    --memory 32

echo "SPAdes assembly complete"
echo "Output files:"
ls -lh $OUTDIR