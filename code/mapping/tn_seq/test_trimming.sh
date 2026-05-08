mkdir -p /proj/uppmax2026-1-61/nobackup/mdha5802/tmp

INDEX=/home/mdha5802/genome_analysis/output/mapping/tnseq_index/E745
FASTQ=/home/mdha5802/1_Zhang_2017/transcriptomics_data/Tn-Seq_BHI/trim_ERR1801012_pass.fastq.gz
TMPDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/tmp

zcat $FASTQ | head -40000 | awk '
NR%4==1 { header=$0; next }
NR%4==2 {
    seq=$0
    if (substr(seq,1,6)=="ATCACG") seq=substr(seq,7)
    pos=index(seq,"TAACAGGTTGGATGATAAGTCCCCGGTCTTC")
    if (pos==0) pos=index(seq,"TAACAGGTTGGATGATAAGT")
    if (pos==0) pos=index(seq,"TAACAGGTTGGAT")
    if (pos==0) pos=index(seq,"TAACAGGTTGG")
    if (pos>0) { read=substr(seq,1,pos-1) } else { read="" }
    next
}
NR%4==3 { qual_header=$0; next }
NR%4==0 {
    qseq=substr($0,7)
    if (length(read)>=14) {
        qual=substr(qseq,1,length(read))
        print header"\n"read"\n"qual_header"\n"qual
    }
}
' > $TMPDIR/test_trimmed.fastq

echo "Reads trimmed: $(wc -l < $TMPDIR/test_trimmed.fastq | awk '{print $1/4}')"

bowtie -v 0 -m 1 --best -p 4 -q -S \
    $INDEX \
    $TMPDIR/test_trimmed.fastq \
    $TMPDIR/test.sam \
    2>&1 | tail -5

samtools view -bS $TMPDIR/test.sam | samtools sort -o $TMPDIR/test.bam
samtools flagstat $TMPDIR/test.bam