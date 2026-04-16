#!/bin/bash
#SBATCH --job-name=rnaseq_mapping
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=rnaseq_mapping_%j.log

module load bowtie2
module load samtools

mkdir -p 1_Zhang_2017/rnaseq/aligned

# Define the index path
INDEX=1_Zhang_2017/rnaseq/genome_index/E745

# Define trimmed data paths
BHI=1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trimmed
SERUM=1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trimmed

# Map BHI samples
for SAMPLE in ERR1797972 ERR1797973 ERR1797974; do
    echo "Mapping BHI sample: $SAMPLE"
    bowtie2 \
        -x $INDEX \
        -1 $BHI/trim_paired_${SAMPLE}_pass_1.fastq.gz \
        -2 $BHI/trim_paired_${SAMPLE}_pass_2.fastq.gz \
        -U $BHI/trim_single_${SAMPLE}_pass_1.fastq.gz,$BHI/trim_single_${SAMPLE}_pass_2.fastq.gz \
        -p 8 \
        --no-unal \
        2> 1_Zhang_2017/rnaseq/aligned/${SAMPLE}_mapping_stats.txt | \
        samtools sort -@ 8 -o 1_Zhang_2017/rnaseq/aligned/${SAMPLE}.bam
    
    samtools index 1_Zhang_2017/rnaseq/aligned/${SAMPLE}.bam
    echo "Done: $SAMPLE"
done

# Map Serum samples
for SAMPLE in ERR1797969 ERR1797970 ERR1797971; do
    echo "Mapping Serum sample: $SAMPLE"
    bowtie2 \
        -x $INDEX \
        -1 $SERUM/trim_paired_${SAMPLE}_pass_1.fastq.gz \
        -2 $SERUM/trim_paired_${SAMPLE}_pass_2.fastq.gz \
        -U $SERUM/trim_single_${SAMPLE}_pass_1.fastq.gz,$SERUM/trim_single_${SAMPLE}_pass_2.fastq.gz \
        -p 8 \
        --no-unal \
        2> 1_Zhang_2017/rnaseq/aligned/${SAMPLE}_mapping_stats.txt | \
        samtools sort -@ 8 -o 1_Zhang_2017/rnaseq/aligned/${SAMPLE}.bam
    
    samtools index 1_Zhang_2017/rnaseq/aligned/${SAMPLE}.bam
    echo "Done: $SAMPLE"
done

echo "All mapping complete"