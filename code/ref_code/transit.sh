#!/bin/bash
#SBATCH --job-name=transit
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=transit_%j.log

module load transit
# or: conda activate transit

mkdir -p 1_Zhang_2017/tnseq/transit_output

# TRANSIT needs the genome in a specific format
# First convert your FASTA to TRANSIT format
python -m pytransit convert fasta \
    1_Zhang_2017/annotation/prokka_output/E745.fasta \
    1_Zhang_2017/tnseq/E745.prot_table

# Run resampling analysis — BHI vs HSerum
python -m pytransit resampling \
    1_Zhang_2017/tnseq/aligned/ERR1801012.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801013.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801014.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801009.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801010.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801011.bam \
    1_Zhang_2017/annotation/prokka_output/E745.gff \
    1_Zhang_2017/tnseq/transit_output/BHI_vs_HSerum.txt \
    -a 0.05

# Run resampling analysis — BHI vs Serum
python -m pytransit resampling \
    1_Zhang_2017/tnseq/aligned/ERR1801012.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801013.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801014.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801006.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801007.bam \
    1_Zhang_2017/tnseq/aligned/ERR1801008.bam \
    1_Zhang_2017/annotation/prokka_output/E745.gff \
    1_Zhang_2017/tnseq/transit_output/BHI_vs_Serum.txt \
    -a 0.05

echo "TRANSIT analysis complete"
