#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J snp_calling
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module load BWA/0.7
module load SAMtools/1.22
module load BCFtools/1.22

# ── Paths ────────────────────────────────────────────────────────
ASSEMBLY=/proj/uppmax2026-1-61/nobackup/mdha5802/assembly/pacbio/E745.contigs.fasta
ILLUMINA=/home/mdha5802/1_Zhang_2017/genomics_data/Illumina
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/output/snp_calling
INDEX=/home/mdha5802/genome_analysis/output/mapping/genome_index/E745

mkdir -p $OUTDIR

# ── Step 1: Index assembly ───────────────────────────────────────
if [ ! -f ${INDEX}.bwt ]; then
    echo "Building BWA index..."
    bwa index -p $INDEX $ASSEMBLY
    echo "Index done"
else
    echo "Index already exists — skipping"
fi

# ── Step 2: Map Illumina reads to assembly ───────────────────────
echo "Mapping Illumina reads..."
bwa mem -t 8 -M \
    $INDEX \
    $ILLUMINA/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
    $ILLUMINA/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
| samtools sort -@ 8 -o $OUTDIR/illumina_mapped.bam

samtools index $OUTDIR/illumina_mapped.bam
samtools flagstat $OUTDIR/illumina_mapped.bam > $OUTDIR/mapping_stats.txt
echo "Mapping stats:"
cat $OUTDIR/mapping_stats.txt

# ── Step 3: Call variants ────────────────────────────────────────
echo "Calling variants..."
bcftools mpileup \
    -f $ASSEMBLY \
    -a AD,DP,SP \
    --min-MQ 20 \
    --min-BQ 20 \
    $OUTDIR/illumina_mapped.bam | \
bcftools call \
    -mv \
    -Oz \
    -o $OUTDIR/variants_raw.vcf.gz

bcftools index $OUTDIR/variants_raw.vcf.gz

# ── Step 4: Filter high quality SNPs ────────────────────────────
echo "Filtering variants..."
bcftools filter \
    -i 'QUAL>30 && FORMAT/DP>10 && INFO/MQ>30' \
    $OUTDIR/variants_raw.vcf.gz \
    -Oz -o $OUTDIR/variants_filtered.vcf.gz

bcftools index $OUTDIR/variants_filtered.vcf.gz

# ── Step 5: Summary statistics ───────────────────────────────────
echo "=== SNP Summary ==="
echo "Raw variants:"
bcftools stats $OUTDIR/variants_raw.vcf.gz | grep "^SN"

echo "Filtered variants:"
bcftools stats $OUTDIR/variants_filtered.vcf.gz | grep "^SN"