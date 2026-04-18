#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J trimmomatic_rnaseq
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load Trimmomatic

# ── Paths ────────────────────────────────────────────────────────
DATADIR=/home/mdha5802/1_Zhang_2017/transcriptomics_data
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/trimmed/transcriptomics
ADAPTERS=/home/mdha5802/genome_analysis/data/trimming/adapters/TruSeq3-PE.fa

mkdir -p /home/mdha5802/genome_analysis/data/trimming/adapters
mkdir -p $OUTDIR/BHI
mkdir -p $OUTDIR/Serum

# ── Download adapter file if not exists ──────────────────────────
if [ ! -f $ADAPTERS ]; then
    echo "Adapter file not found — downloading..."
    wget https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa \
        -O $ADAPTERS
    echo "Adapter file downloaded: $ADAPTERS"
else
    echo "Adapter file found — skipping download"
fi

# ── Quick test on 10,000 reads before full run ───────────────────
# echo "=== Running quick test ==="
# TESTDIR=$OUTDIR/test/trimming
# mkdir -p $TESTDIR

# zcat $DATADIR/RNA-Seq_BH/raw/ERR1797972_1.fastq.gz | head -40000 | gzip > $TESTDIR/test_R1.fastq.gz
# zcat $DATADIR/RNA-Seq_BH/raw/ERR1797972_2.fastq.gz | head -40000 | gzip > $TESTDIR/test_R2.fastq.gz

# trimmomatic PE -threads 4 \
#     $TESTDIR/test_R1.fastq.gz \
#     $TESTDIR/test_R2.fastq.gz \
#     $TESTDIR/out_paired_1.fastq.gz \
#     $TESTDIR/out_single_1.fastq.gz \
#     $TESTDIR/out_paired_2.fastq.gz \
#     $TESTDIR/out_single_2.fastq.gz \
#     ILLUMINACLIP:$ADAPTERS:2:30:10 \
#     LEADING:3 TRAILING:3 \
#     SLIDINGWINDOW:4:15 \
#     MINLEN:36

# if [ $? -eq 0 ]; then
#     echo "Test passed — proceeding with full trimming"
#     rm -rf $TESTDIR
# else
#     echo "ERROR: Test failed — check parameters before full run"
#     exit 1
# fi

# ── Trim BHI samples ─────────────────────────────────────────────
echo "=== Trimming BHI samples ==="

for SAMPLE in ERR1797972 ERR1797973 ERR1797974; do
    echo "Trimming BHI: $SAMPLE"
    trimmomatic PE -threads 4 \
        $DATADIR/RNA-Seq_BH/raw/${SAMPLE}_1.fastq.gz \
        $DATADIR/RNA-Seq_BH/raw/${SAMPLE}_2.fastq.gz \
        $OUTDIR/BHI/trim_paired_${SAMPLE}_pass_1.fastq.gz \
        $OUTDIR/BHI/trim_single_${SAMPLE}_pass_1.fastq.gz \
        $OUTDIR/BHI/trim_paired_${SAMPLE}_pass_2.fastq.gz \
        $OUTDIR/BHI/trim_single_${SAMPLE}_pass_2.fastq.gz \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
    echo "Done: $SAMPLE"
done

# ── Trim Serum samples ───────────────────────────────────────────
echo "=== Trimming Serum samples ==="

for SAMPLE in ERR1797969 ERR1797970 ERR1797971; do
    echo "Trimming Serum: $SAMPLE"
    trimmomatic PE -threads 4 \
        $DATADIR/RNA-Seq_Serum/raw/${SAMPLE}_1.fastq.gz \
        $DATADIR/RNA-Seq_Serum/raw/${SAMPLE}_2.fastq.gz \
        $OUTDIR/Serum/trim_paired_${SAMPLE}_pass_1.fastq.gz \
        $OUTDIR/Serum/trim_single_${SAMPLE}_pass_1.fastq.gz \
        $OUTDIR/Serum/trim_paired_${SAMPLE}_pass_2.fastq.gz \
        $OUTDIR/Serum/trim_single_${SAMPLE}_pass_2.fastq.gz \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
    echo "Done: $SAMPLE"
done

echo "=== All trimming complete ==="
echo "BHI output:"
ls -lh $OUTDIR/BHI
echo "Serum output:"
ls -lh $OUTDIR/Serum