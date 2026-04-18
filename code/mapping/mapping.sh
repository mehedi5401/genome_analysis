#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 10:00:00
#SBATCH -J bwa_mapping
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load BWA/0.7
module load SAMtools/1.22

# ── Paths ────────────────────────────────────────────────────────
GENOME=/home/mdha5802/genome_analysis/output/annotation/E745.fna
INDEX=/home/mdha5802/genome_analysis/output/mapping/genome_index/E745
TRIMMED=/proj/uppmax2026-1-61/nobackup/mdha5802/trimmed/transcriptomics
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/rna_seq
TMPDIR=$OUTDIR/TMP

# Set TEST=true to run on 10,000 reads only — set to false for full run
TEST=false

mkdir -p /home/mdha5802/genome_analysis/output/mapping/genome_index
mkdir -p $OUTDIR/BHI
mkdir -p $OUTDIR/Serum

# ── Step 1: Build genome index (once only) ───────────────────────
if [ ! -f ${INDEX}.bwt ]; then
    echo "Building BWA index..."
    bwa index -p $INDEX $GENOME
    echo "Index done"
else
    echo "Index already exists — skipping"
fi

# ── Step 2: Map BHI samples ──────────────────────────────────────
echo "=== Mapping BHI samples ==="

for SAMPLE in ERR1797972 ERR1797973 ERR1797974; do
    echo "Mapping BHI: $SAMPLE"

    if [ "$TEST" = true ]; then
        mkdir -p $TMPDIR
        zcat $TRIMMED/BHI/trim_paired_${SAMPLE}_pass_1.fastq.gz | head -40000 | gzip > $TMPDIR/${SAMPLE}_R1.fastq.gz
        zcat $TRIMMED/BHI/trim_paired_${SAMPLE}_pass_2.fastq.gz | head -40000 | gzip > $TMPDIR/${SAMPLE}_R2.fastq.gz
        zcat $TRIMMED/BHI/trim_single_${SAMPLE}_pass_1.fastq.gz | head -40000 | gzip > $TMPDIR/${SAMPLE}_S1.fastq.gz
        zcat $TRIMMED/BHI/trim_single_${SAMPLE}_pass_2.fastq.gz | head -40000 | gzip > $TMPDIR/${SAMPLE}_S2.fastq.gz
        R1=$TMPDIR/${SAMPLE}_R1.fastq.gz
        R2=$TMPDIR/${SAMPLE}_R2.fastq.gz
        S1=$TMPDIR/${SAMPLE}_S1.fastq.gz
        S2=$TMPDIR/${SAMPLE}_S2.fastq.gz
    else
        R1=$TRIMMED/BHI/trim_paired_${SAMPLE}_pass_1.fastq.gz
        R2=$TRIMMED/BHI/trim_paired_${SAMPLE}_pass_2.fastq.gz
        S1=$TRIMMED/BHI/trim_single_${SAMPLE}_pass_1.fastq.gz
        S2=$TRIMMED/BHI/trim_single_${SAMPLE}_pass_2.fastq.gz
    fi

    bwa mem -t 4 -M \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:BHI\tPL:ILLUMINA" \
        $INDEX $R1 $R2 \
    | samtools sort -@ 4 -o $OUTDIR/BHI/${SAMPLE}_paired.bam

    bwa mem -t 4 -M \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:BHI\tPL:ILLUMINA" \
        $INDEX $S1 \
    | samtools sort -@ 4 -o $OUTDIR/BHI/${SAMPLE}_single_1.bam

    bwa mem -t 4 -M \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:BHI\tPL:ILLUMINA" \
        $INDEX $S2 \
    | samtools sort -@ 4 -o $OUTDIR/BHI/${SAMPLE}_single_2.bam

    samtools merge -f \
        $OUTDIR/BHI/${SAMPLE}.bam \
        $OUTDIR/BHI/${SAMPLE}_paired.bam \
        $OUTDIR/BHI/${SAMPLE}_single_1.bam \
        $OUTDIR/BHI/${SAMPLE}_single_2.bam

    samtools index $OUTDIR/BHI/${SAMPLE}.bam
    samtools flagstat $OUTDIR/BHI/${SAMPLE}.bam \
        > $OUTDIR/BHI/${SAMPLE}_flagstat.txt

    rm $OUTDIR/BHI/${SAMPLE}_paired.bam \
       $OUTDIR/BHI/${SAMPLE}_single_1.bam \
       $OUTDIR/BHI/${SAMPLE}_single_2.bam

    echo "Done BHI $SAMPLE: $(grep 'mapped (' $OUTDIR/BHI/${SAMPLE}_flagstat.txt)"
done

# ── Step 3: Map Serum samples ────────────────────────────────────
echo "=== Mapping Serum samples ==="

for SAMPLE in ERR1797969 ERR1797970 ERR1797971; do
    echo "Mapping Serum: $SAMPLE"

    if [ "$TEST" = true ]; then
        mkdir -p $TMPDIR
        zcat $TRIMMED/Serum/trim_paired_${SAMPLE}_pass_1.fastq.gz | head -40000 | gzip > $TMPDIR/${SAMPLE}_R1.fastq.gz
        zcat $TRIMMED/Serum/trim_paired_${SAMPLE}_pass_2.fastq.gz | head -40000 | gzip > $TMPDIR/${SAMPLE}_R2.fastq.gz
        zcat $TRIMMED/Serum/trim_single_${SAMPLE}_pass_1.fastq.gz | head -40000 | gzip > $TMPDIR/${SAMPLE}_S1.fastq.gz
        zcat $TRIMMED/Serum/trim_single_${SAMPLE}_pass_2.fastq.gz | head -40000 | gzip > $TMPDIR/${SAMPLE}_S2.fastq.gz
        R1=$TMPDIR/${SAMPLE}_R1.fastq.gz
        R2=$TMPDIR/${SAMPLE}_R2.fastq.gz
        S1=$TMPDIR/${SAMPLE}_S1.fastq.gz
        S2=$TMPDIR/${SAMPLE}_S2.fastq.gz
    else
        R1=$TRIMMED/Serum/trim_paired_${SAMPLE}_pass_1.fastq.gz
        R2=$TRIMMED/Serum/trim_paired_${SAMPLE}_pass_2.fastq.gz
        S1=$TRIMMED/Serum/trim_single_${SAMPLE}_pass_1.fastq.gz
        S2=$TRIMMED/Serum/trim_single_${SAMPLE}_pass_2.fastq.gz
    fi

    bwa mem -t 4 -M \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:Serum\tPL:ILLUMINA" \
        $INDEX $R1 $R2 \
    | samtools sort -@ 4 -o $OUTDIR/Serum/${SAMPLE}_paired.bam

    bwa mem -t 4 -M \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:Serum\tPL:ILLUMINA" \
        $INDEX $S1 \
    | samtools sort -@ 4 -o $OUTDIR/Serum/${SAMPLE}_single_1.bam

    bwa mem -t 4 -M \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:Serum\tPL:ILLUMINA" \
        $INDEX $S2 \
    | samtools sort -@ 4 -o $OUTDIR/Serum/${SAMPLE}_single_2.bam

    samtools merge -f \
        $OUTDIR/Serum/${SAMPLE}.bam \
        $OUTDIR/Serum/${SAMPLE}_paired.bam \
        $OUTDIR/Serum/${SAMPLE}_single_1.bam \
        $OUTDIR/Serum/${SAMPLE}_single_2.bam

    samtools index $OUTDIR/Serum/${SAMPLE}.bam
    samtools flagstat $OUTDIR/Serum/${SAMPLE}.bam \
        > $OUTDIR/Serum/${SAMPLE}_flagstat.txt

    rm $OUTDIR/Serum/${SAMPLE}_paired.bam \
       $OUTDIR/Serum/${SAMPLE}_single_1.bam \
       $OUTDIR/Serum/${SAMPLE}_single_2.bam

    echo "Done Serum $SAMPLE: $(grep 'mapped (' $OUTDIR/Serum/${SAMPLE}_flagstat.txt)"
done

if [ "$TEST" = true ]; then rm -rf $TMPDIR; fi

echo "=== All mapping complete ==="
echo "BAM files:"
ls -lh $OUTDIR/BHI/*.bam $OUTDIR/Serum/*.bam