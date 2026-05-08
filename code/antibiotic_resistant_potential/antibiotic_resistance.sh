#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH -J resfinder
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# ── Load conda cleanly ───────────────────────────────────────────
module purge
module load Miniconda3/25.7.0-2
conda activate /proj/uppmax2026-1-61/nobackup/mdha5802/envs/resfinder

# ── Paths ────────────────────────────────────────────────────────
ASSEMBLY=/proj/uppmax2026-1-61/nobackup/mdha5802/assembly/pacbio/E745.contigs.fasta
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/output/antibiotic_resistance
DB_RES=/proj/uppmax2026-1-61/nobackup/mdha5802/antibiotic_resistance/db/resfinder_db/
DB_POINT=/proj/uppmax2026-1-61/nobackup/mdha5802/antibiotic_resistance/db/pointfinder_db/

mkdir -p $OUTDIR

# ── Run ResFinder ────────────────────────────────────────────────
echo "Running ResFinder..."

python -m resfinder \
    -ifa $ASSEMBLY \
    -o $OUTDIR \
    -s "Enterococcus faecium" \
    -l 0.6 \
    -t 0.8 \
    --acquired \
    --point \
    --db_path_res $DB_RES \
    --db_path_point $DB_POINT

echo "ResFinder complete"
echo "Output files:"
ls -lh $OUTDIR
echo ""
echo "=== Resistance summary ==="
cat $OUTDIR/ResFinder_results_tab.txt