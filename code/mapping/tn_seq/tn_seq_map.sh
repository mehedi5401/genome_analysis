#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 08:00:00
#SBATCH -J tnseq_mapping
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module purge
module load Bowtie/1.3
module load SAMtools/1.22

# ── Paths ─────────────────────────────────────────────────────────
GENOME=/home/mdha5802/genome_analysis/output/annotation/E745.fna
INDEX=/home/mdha5802/genome_analysis/output/mapping/tnseq_index/E745
DATADIR=/home/mdha5802/1_Zhang_2017/transcriptomics_data
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped
TRIMDIR=$OUTDIR/trimmed
TMPDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/tmp

mkdir -p $OUTDIR/BHI $OUTDIR/HSerum $OUTDIR/Serum $TRIMDIR $TMPDIR
mkdir -p $(dirname $INDEX)

# ── Build bowtie index ────────────────────────────────────────────
if [ ! -f ${INDEX}.1.ebwt ]; then
    echo "Building bowtie index..."
    bowtie-build $GENOME $INDEX \
        || { echo "ERROR: bowtie-build failed"; exit 1; }
    echo "Index done"
else
    echo "Index exists — skipping"
fi

# ── Trim function (barcode-aware) ─────────────────────────────────
trim_reads() {
    local FASTQ=$1
    local TRIMMED=$2
    local BARCODE=$3
    local BARCLEN=${#BARCODE}

    zcat $FASTQ | awk -v bc="$BARCODE" -v bclen="$BARCLEN" '
    NR%4==1 { header=$0; next }
    NR%4==2 {
        seq=$0
        # strip sample-specific barcode from start
        if (substr(seq,1,bclen)==bc) seq=substr(seq,bclen+1)
        # find transposon junction
        pos=index(seq,"TAACAGGTTGGATGATAAGTCCCCGGTCTTC")
        if (pos==0) pos=index(seq,"TAACAGGTTGGATGATAAGT")
        if (pos==0) pos=index(seq,"TAACAGGTTGGAT")
        if (pos==0) pos=index(seq,"TAACAGGTTGG")
        if (pos>0) { read=substr(seq,1,pos-1) } else { read="" }
        next
    }
    NR%4==3 { qual_header=$0; next }
    NR%4==0 {
        qseq=substr($0,bclen+1)
        if (length(read)>=14) {
            qual=substr(qseq,1,length(read))
            print header"\n"read"\n"qual_header"\n"qual
        }
    }
    ' > $TRIMMED
    echo "  Reads after trimming: $(wc -l < $TRIMMED | awk '{print $1/4}')"
}

# ── Map function ──────────────────────────────────────────────────
process_sample() {
    local SAMPLE=$1
    local FASTQ=$2
    local OUTBAM=$3
    local BARCODE=$4

    echo "=== Processing $SAMPLE (barcode: $BARCODE) ==="

    # Step 1: Trim
    TRIMMED=$TRIMDIR/${SAMPLE}_trimmed.fastq
    echo "  Trimming..."
    trim_reads $FASTQ $TRIMMED $BARCODE

    # Step 2: Map with bowtie
    SAM=$TMPDIR/${SAMPLE}.sam
    echo "  Mapping..."
    bowtie \
        -v 0 \
        -m 1 \
        --best \
        -p 4 \
        -q \
        -S \
        $INDEX \
        $TRIMMED \
        $SAM \
        2> $OUTDIR/${SAMPLE}_bowtie.log \
        || { echo "ERROR: bowtie failed for $SAMPLE"; exit 1; }
    cat $OUTDIR/${SAMPLE}_bowtie.log

    # Step 3: SAM to BAM
    echo "  Converting to BAM..."
    samtools view -bS $SAM \
        | samtools sort -@ 4 -o $OUTBAM \
        || { echo "ERROR: samtools failed for $SAMPLE"; exit 1; }
    samtools index $OUTBAM

    # cleanup
    rm $SAM $TRIMMED

    # Step 4: Stats
    echo "  Stats for $SAMPLE:"
    samtools flagstat $OUTBAM | tee ${OUTBAM%.bam}_flagstat.txt

    echo "=== Done: $SAMPLE ==="
}

# ── BHI ───────────────────────────────────────────────────────────
echo "=== Mapping BHI samples ==="
process_sample BHI_ERR1801012 $DATADIR/Tn-Seq_BHI/trim_ERR1801012_pass.fastq.gz $OUTDIR/BHI/ERR1801012.bam ATCACG
process_sample BHI_ERR1801013 $DATADIR/Tn-Seq_BHI/trim_ERR1801013_pass.fastq.gz $OUTDIR/BHI/ERR1801013.bam ACAGTG
process_sample BHI_ERR1801014 $DATADIR/Tn-Seq_BHI/trim_ERR1801014_pass.fastq.gz $OUTDIR/BHI/ERR1801014.bam GCCAAT

# ── HSerum ────────────────────────────────────────────────────────
echo "=== Mapping HSerum samples ==="
process_sample HSerum_ERR1801009 $DATADIR/Tn-Seq_HSerum/trim_ERR1801009_pass.fastq.gz $OUTDIR/HSerum/ERR1801009.bam CGATGT
process_sample HSerum_ERR1801010 $DATADIR/Tn-Seq_HSerum/trim_ERR1801010_pass.fastq.gz $OUTDIR/HSerum/ERR1801010.bam TTAGGC
process_sample HSerum_ERR1801011 $DATADIR/Tn-Seq_HSerum/trim_ERR1801011_pass.fastq.gz $OUTDIR/HSerum/ERR1801011.bam CAGATC

# ── Serum ─────────────────────────────────────────────────────────
echo "=== Mapping Serum samples ==="
process_sample Serum_ERR1801006 $DATADIR/Tn-Seq_Serum/trim_ERR1801006_pass.fastq.gz $OUTDIR/Serum/ERR1801006.bam ACTTGA
process_sample Serum_ERR1801007 $DATADIR/Tn-Seq_Serum/trim_ERR1801007_pass.fastq.gz $OUTDIR/Serum/ERR1801007.bam GATCAG
process_sample Serum_ERR1801008 $DATADIR/Tn-Seq_Serum/trim_ERR1801008_pass.fastq.gz $OUTDIR/Serum/ERR1801008.bam TAGCTT

echo "=== All mapping complete ==="
echo "BAM files:"
ls -lh $OUTDIR/BHI/ $OUTDIR/HSerum/ $OUTDIR/Serum/

echo "=== Mapping rate summary ==="
for BAM in $OUTDIR/BHI/*.bam $OUTDIR/HSerum/*.bam $OUTDIR/Serum/*.bam; do
    echo "--- $(basename $BAM) ---"
    samtools flagstat $BAM | grep "mapped ("
done