#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=featurecounts_%j.log

module load subread

mkdir -p 1_Zhang_2017/rnaseq/counts

featureCounts \
    -a 1_Zhang_2017/annotation/prokka_output/E745.gff \
    -o 1_Zhang_2017/rnaseq/counts/all_samples_counts.txt \
    -t CDS \
    -g ID \
    -p \
    -T 8 \
    1_Zhang_2017/rnaseq/aligned/ERR1797972.bam \
    1_Zhang_2017/rnaseq/aligned/ERR1797973.bam \
    1_Zhang_2017/rnaseq/aligned/ERR1797974.bam \
    1_Zhang_2017/rnaseq/aligned/ERR1797969.bam \
    1_Zhang_2017/rnaseq/aligned/ERR1797970.bam \
    1_Zhang_2017/rnaseq/aligned/ERR1797971.bam

echo "Counting complete"