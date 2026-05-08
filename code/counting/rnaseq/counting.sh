#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 02:00:00
#SBATCH -J htseq_count
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

module load HTSeq/2.1
module load SAMtools/1.22

# ── Paths ────────────────────────────────────────────────────────
GFF=/home/mdha5802/genome_analysis/output/annotation/E745_clean.gff
BAMDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/rna_seq
OUTDIR=/home/mdha5802/genome_analysis/output/counting/rnaseq

mkdir -p $OUTDIR

# Set TEST=true to run on first BAM only for quick check
TEST=false

# ── Count BHI samples ────────────────────────────────────────────
echo "=== Counting BHI samples ==="

for SAMPLE in ERR1797972 ERR1797973 ERR1797974; do
    echo "Counting BHI: $SAMPLE"

    # Filter to properly paired primary reads only
    samtools view -b -f 0x2 -F 0x100 \
        $BAMDIR/BHI/${SAMPLE}.bam \
        > $BAMDIR/BHI/${SAMPLE}_clean.bam
    samtools index $BAMDIR/BHI/${SAMPLE}_clean.bam

    htseq-count \
        --format bam \
        --order pos \
        --stranded no \
        --type CDS \
        --idattr ID \
        --quiet \
        $BAMDIR/BHI/${SAMPLE}_clean.bam \
        $GFF \
        > $OUTDIR/${SAMPLE}_counts.txt \
        2> $OUTDIR/${SAMPLE}_counts.log

    rm $BAMDIR/BHI/${SAMPLE}_clean.bam $BAMDIR/BHI/${SAMPLE}_clean.bam.bai

    if [ -s $OUTDIR/${SAMPLE}_counts.txt ]; then
        echo "Done: $SAMPLE — $(grep -v '^__' $OUTDIR/${SAMPLE}_counts.txt | wc -l) genes counted"
    else
        echo "ERROR: $SAMPLE — check $OUTDIR/${SAMPLE}_counts.log"
        cat $OUTDIR/${SAMPLE}_counts.log
    fi

    if [ "$TEST" = true ]; then
        echo "Test mode — stopping after first sample"
        break
    fi
done

# ── Count Serum samples ──────────────────────────────────────────
echo "=== Counting Serum samples ==="

for SAMPLE in ERR1797969 ERR1797970 ERR1797971; do
    echo "Counting Serum: $SAMPLE"

    # Filter to properly paired primary reads only
    samtools view -b -f 0x2 -F 0x100 \
        $BAMDIR/Serum/${SAMPLE}.bam \
        > $BAMDIR/Serum/${SAMPLE}_clean.bam
    samtools index $BAMDIR/Serum/${SAMPLE}_clean.bam

    htseq-count \
        --format bam \
        --order pos \
        --stranded no \
        --type CDS \
        --idattr ID \
        --quiet \
        $BAMDIR/Serum/${SAMPLE}_clean.bam \
        $GFF \
        > $OUTDIR/${SAMPLE}_counts.txt \
        2> $OUTDIR/${SAMPLE}_counts.log

    rm $BAMDIR/Serum/${SAMPLE}_clean.bam $BAMDIR/Serum/${SAMPLE}_clean.bam.bai

    if [ -s $OUTDIR/${SAMPLE}_counts.txt ]; then
        echo "Done: $SAMPLE — $(grep -v '^__' $OUTDIR/${SAMPLE}_counts.txt | wc -l) genes counted"
    else
        echo "ERROR: $SAMPLE — check $OUTDIR/${SAMPLE}_counts.log"
        cat $OUTDIR/${SAMPLE}_counts.log
    fi

    if [ "$TEST" = true ]; then
        echo "Test mode — stopping after first sample"
        break
    fi
done

echo "=== All counting complete ==="
echo "Output files:"
ls -lh $OUTDIR/

# ── Quick sanity check ───────────────────────────────────────────
echo "=== Summary statistics ==="
for FILE in $OUTDIR/*_counts.txt; do
    SAMPLE=$(basename $FILE _counts.txt)
    TOTAL=$(grep -v '^__' $FILE | awk '{sum+=$2} END {print sum}')
    NO_FEAT=$(grep '__no_feature' $FILE | awk '{print $2}')
    echo "$SAMPLE — total counts: $TOTAL — no_feature: $NO_FEAT"
done