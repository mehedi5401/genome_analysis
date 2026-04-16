#!/bin/bash
#SBATCH --job-name=tnseq_mapping
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=tnseq_mapping_%j.log

module load bowtie2
module load samtools

mkdir -p 1_Zhang_2017/tnseq/aligned

INDEX=1_Zhang_2017/rnaseq/genome_index/E745

# Map BHI samples
for SAMPLE in ERR1801012 ERR1801013 ERR1801014; do
    echo "Mapping Tn-seq BHI: $SAMPLE"
    bowtie2 \
        -x $INDEX \
        -U 1_Zhang_2017/transcriptomics_data/Tn-Seq_BHI/trim_${SAMPLE}_pass.fastq.gz \
        -p 8 \
        --no-unal \
        --very-sensitive \
        2> 1_Zhang_2017/tnseq/aligned/${SAMPLE}_stats.txt | \
        samtools sort -@ 8 -o 1_Zhang_2017/tnseq/aligned/${SAMPLE}.bam
    samtools index 1_Zhang_2017/tnseq/aligned/${SAMPLE}.bam
done

# Map HSerum samples
for SAMPLE in ERR1801009 ERR1801010 ERR1801011; do
    echo "Mapping Tn-seq HSerum: $SAMPLE"
    bowtie2 \
        -x $INDEX \
        -U 1_Zhang_2017/transcriptomics_data/Tn-Seq_HSerum/trim_${SAMPLE}_pass.fastq.gz \
        -p 8 \
        --no-unal \
        --very-sensitive \
        2> 1_Zhang_2017/tnseq/aligned/${SAMPLE}_stats.txt | \
        samtools sort -@ 8 -o 1_Zhang_2017/tnseq/aligned/${SAMPLE}.bam
    samtools index 1_Zhang_2017/tnseq/aligned/${SAMPLE}.bam
done

# Map Serum samples
for SAMPLE in ERR1801006 ERR1801007 ERR1801008; do
    echo "Mapping Tn-seq Serum: $SAMPLE"
    bowtie2 \
        -x $INDEX \
        -U 1_Zhang_2017/transcriptomics_data/Tn-Seq_Serum/trim_${SAMPLE}_pass.fastq.gz \
        -p 8 \
        --no-unal \
        --very-sensitive \
        2> 1_Zhang_2017/tnseq/aligned/${SAMPLE}_stats.txt | \
        samtools sort -@ 8 -o 1_Zhang_2017/tnseq/aligned/${SAMPLE}.bam
    samtools index 1_Zhang_2017/tnseq/aligned/${SAMPLE}.bam
done

echo "All Tn-seq mapping complete"