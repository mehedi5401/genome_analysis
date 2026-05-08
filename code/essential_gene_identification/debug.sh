PYTHON=/proj/uppmax2026-1-61/nobackup/mdha5802/envs/resfinder/bin/python
IDS=/home/mdha5802/genome_analysis/output/essential_gene_identification/E745_TA_sites.txt

$PYTHON - <<EOF
import pysam

ta_sites = set()
with open("$IDS") as f:
    for line in f:
        chrom, pos = line.strip().split("\t")
        ta_sites.add((chrom, int(pos)))

samples = {
    "BHI_ERR1801012":    "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/BHI/ERR1801012.bam",
    "BHI_ERR1801013":    "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/BHI/ERR1801013.bam",
    "BHI_ERR1801014":    "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/BHI/ERR1801014.bam",
    "HSerum_ERR1801009": "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/HSerum/ERR1801009.bam",
    "HSerum_ERR1801010": "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/HSerum/ERR1801010.bam",
    "HSerum_ERR1801011": "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/HSerum/ERR1801011.bam",
    "Serum_ERR1801006":  "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/Serum/ERR1801006.bam",
    "Serum_ERR1801007":  "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/Serum/ERR1801007.bam",
    "Serum_ERR1801008":  "/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped/Serum/ERR1801008.bam",
}

for sample, bam_path in samples.items():
    counts = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam:
        if read.is_unmapped:
            continue
        chrom = read.reference_name
        pos = read.reference_end if read.is_reverse else read.reference_start + 1
        if (chrom, pos) in ta_sites:
            counts[(chrom, pos)] = counts.get((chrom, pos), 0) + 1
    bam.close()
    unique = len(counts)
    total = sum(counts.values())
    avg = total/unique if unique > 0 else 0
    print(f"{sample}: {unique} unique sites, {total} total counts, avg depth {avg:.0f}x")
EOF