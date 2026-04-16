#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
```
#temporary chagne to run in login node as uppmax was not available.
# SBATCH -p pelle
# SBATCH -c 2
# SBATCH -t 00:30:00
# SBATCH -J pacbio_assembly_evaluation
# SBATCH --mail-type=ALL
# SBATCH --output=%x.%j.out
```

# Load Prokka module (adjust to your HPC environment)
module load prokka/1.14

# ── Input / output ──────────────────────────────────────────
GENOME=/proj/uppmax2026-1-61/nobackup/mdha5802/assembly/pacbio/E745.contigs.fasta
OUTDIR=/home/mdha5802/genome_analysis/output/annotation
# ── Run Prokka ───────────────────────────────────────────────
prokka $GENOME \
    --outdir  $OUTDIR \
    --prefix E745 \
    --genus Enterococcus \
    --species faecium \
    --cpus 2 \
    --mincontiglen 500 \
    --kingdom Bacteria \
    --force
    

echo "Prokka complete"
echo "Output files:"
ls -lh $OUTDIR

