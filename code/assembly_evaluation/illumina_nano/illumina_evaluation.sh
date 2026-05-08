#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
# SBATCH -p pelle
# SBATCH -c 2
# SBATCH -t 01:00:00
# SBATCH -J illumina_nano_assembly_evaluation
# SBATCH --mail-type=ALL
# SBATCH --output=%x.%j.out

module load QUAST/5.3
module load BUSCO/5.8


# ── Paths ────────────────────────────────────────────────────────
ASSEMBLY=/proj/uppmax2026-1-61/nobackup/mdha5802/assembly/illumina_nano/contigs.fasta
REFERENCE=/home/mdha5802/genome_analysis/data/assembly_evaluation/mummerplot_ref/E_faecium_reference.fasta
EVALDIR=/home/mdha5802/genome_analysis/output/assembly_evaluation/illumina_nano

mkdir -p $EVALDIR/quast
mkdir -p $EVALDIR/busco
mkdir -p $EVALDIR/mummer

# ── Download reference if not exists ─────────────────────────────
if [ ! -f $REFERENCE ]; then
    echo "Reference not found — downloading from NCBI..."
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP014529,CP014530,CP014531,CP014532,CP014533,CP014534,CP014535&rettype=fasta&retmode=text" \
        -O $REFERENCE
    echo "Reference downloaded: $(grep -c '^>' $REFERENCE) sequences"
else
    echo "Reference already exists — skipping download"
fi

#── QUAST ─────────────────────────────────────────────────────────
echo "Running QUAST..."
quast \
    $ASSEMBLY \
    --reference $REFERENCE \
    --output-dir $EVALDIR/quast \
    --threads 2 \
    --min-contig 500
echo "QUAST done"

#── BUSCO ─────────────────────────────────────────────────────────
echo "Running BUSCO..."
busco \
    -i $ASSEMBLY \
    -o busco_results \
    --out_path $EVALDIR/busco \
    -m genome \
    -l lactobacillales_odb10 \
    --cpu 2 \
    -f
echo "BUSCO done"



# ── MUMmer + MUMmerplot ───────────────────────────────────────────
module purge
module load MUMmer/4.0

echo "Running MUMmer alignment..."
nucmer \
    --prefix $EVALDIR/mummer/E745_vs_ref \
    --threads 2 \
    $REFERENCE \
    $ASSEMBLY

echo "Filtering alignments..."
delta-filter \
    -r -q \
    $EVALDIR/mummer/E745_vs_ref.delta \
    > $EVALDIR/mummer/E745_vs_ref.filtered.delta

echo "Drawing dot plot..."
mummerplot \
    $EVALDIR/mummer/E745_vs_ref.filtered.delta \
    --png \
    --large \
    -p $EVALDIR/mummer/E745_vs_ref_dotplot
echo "MUMmer done"

echo "All evaluation complete"
echo "Results in: $EVALDIR"