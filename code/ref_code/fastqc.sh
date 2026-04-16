#!/bin/bash
# FastQC on all raw reads
fastqc \
    genomics_data/Illumina/*.fq.gz \
    transcriptomics_data/RNA-Seq_BH/raw/*.fastq.gz \
    transcriptomics_data/RNA-Seq_Serum/raw/*.fastq.gz \
    --outdir QC/fastqc/raw/ \
    --threads 4

# Summarise all reports
multiqc QC/fastqc/raw/ --outdir QC/multiqc/raw/