#!/bin/bash
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 8
#SBATCH -t 10:00:00
#SBATCH --mem=32G
#SBATCH -J eggnog_E745
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-type=ALL

module load eggnog-mapper/2.1

INDIR=/home/mdha5802/genome_analysis/output/annotation
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/output/refined_annotation
DB=/proj/uppmax2026-1-61/nobackup/work/qich5654/eggnog_db/
mkdir -p $OUTDIR

emapper.py \
  -i $INDIR/E745.faa \
  --itype proteins \
  -o E745 \
  --output_dir $OUTDIR \
  --cpu 2 \
  --data_dir $DB \
  --tax_scope Firmicutes \
  --override