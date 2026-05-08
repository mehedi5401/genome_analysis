#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 08:00:00
#SBATCH -J transit
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module purge
module load Miniconda3/25.7.0-2
conda init bash
source ~/.bashrc
conda activate /proj/uppmax2026-1-61/nobackup/mdha5802/envs/resfinder

PYTHON=/proj/uppmax2026-1-61/nobackup/mdha5802/envs/resfinder/bin/python

# ── Verify pytransit ──────────────────────────────────────────────
$PYTHON -m pytransit -version || { echo "pytransit not found, exiting"; exit 1; }

# ── Paths ─────────────────────────────────────────────────────────
FASTA=/home/mdha5802/genome_analysis/output/annotation/E745.fna
GFF=/home/mdha5802/genome_analysis/output/annotation/E745_clean.gff
BAMDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/mapping/tnseq_remapped
OUTDIR=/home/mdha5802/genome_analysis/output/essential_gene_identification
WIGDIR=$OUTDIR/wigs
FILTERED=$WIGDIR/filtered
PROT_TABLE=$OUTDIR/E745.prot_table  # ← add this line

mkdir -p $OUTDIR $WIGDIR $FILTERED

# # ── Step 1: Convert GFF to prot_table ────────────────────────────
# echo "Converting GFF to prot_table..."
# PROT_TABLE=$OUTDIR/E745.prot_table
# $PYTHON -m pytransit convert gff_to_prot \
#     $GFF \
#     $PROT_TABLE \
#     || { echo "ERROR: gff_to_prot failed"; exit 1; }
# echo "prot_table done: $(wc -l < $PROT_TABLE) genes"

# # ── Step 2: Generate TA site ids from FASTA ───────────────────────
# echo "Generating TA sites from FASTA..."
# IDS_FILE=$OUTDIR/E745_TA_sites.txt
# $PYTHON - <<EOF
# import re

# fasta = {}
# current = None
# with open("$FASTA") as f:
#     for line in f:
#         line = line.strip()
#         if line.startswith(">"):
#             current = line[1:].split()[0]
#             fasta[current] = []
#         else:
#             fasta[current].append(line)

# with open("$IDS_FILE", "w") as out:
#     for chrom, seqs in fasta.items():
#         seq = "".join(seqs).upper()
#         for m in re.finditer("TA", seq):
#             out.write(f"{chrom}\t{m.start()+1}\n")

# print(f"TA sites written: $(wc -l < $IDS_FILE)")
# EOF
# echo "TA sites found: $(wc -l < $IDS_FILE)"

# # ── Step 3: BAM to WIG using pysam ───────────────────────────────
# echo "Converting BAM files to WIG..."
# $PYTHON - <<EOF
# import pysam

# IDS_FILE = "$IDS_FILE"
# BAMDIR   = "$BAMDIR"
# WIGDIR   = "$WIGDIR"

# # load TA sites
# print("Loading TA sites...")
# ta_sites = {}  # chrom -> sorted list of positions
# with open(IDS_FILE) as f:
#     for line in f:
#         chrom, pos = line.strip().split("\t")
#         pos = int(pos)
#         if chrom not in ta_sites:
#             ta_sites[chrom] = []
#         ta_sites[chrom].append(pos)

# # sort positions per chrom
# for chrom in ta_sites:
#     ta_sites[chrom].sort()

# total_sites = sum(len(v) for v in ta_sites.values())
# print(f"Loaded {total_sites} TA sites across {len(ta_sites)} contigs")

# # define samples
# samples = {
#     "BHI_ERR1801012":    f"{BAMDIR}/BHI/ERR1801012.bam",
#     "BHI_ERR1801013":    f"{BAMDIR}/BHI/ERR1801013.bam",
#     "BHI_ERR1801014":    f"{BAMDIR}/BHI/ERR1801014.bam",
#     "HSerum_ERR1801009": f"{BAMDIR}/HSerum/ERR1801009.bam",
#     "HSerum_ERR1801010": f"{BAMDIR}/HSerum/ERR1801010.bam",
#     "HSerum_ERR1801011": f"{BAMDIR}/HSerum/ERR1801011.bam",
#     "Serum_ERR1801006":  f"{BAMDIR}/Serum/ERR1801006.bam",
#     "Serum_ERR1801007":  f"{BAMDIR}/Serum/ERR1801007.bam",
#     "Serum_ERR1801008":  f"{BAMDIR}/Serum/ERR1801008.bam",
# }

# for sample, bam_path in samples.items():
#     print(f"Processing {sample}...")
#     counts = {}

#     bam = pysam.AlignmentFile(bam_path, "rb")
#     for read in bam:
#         if read.is_unmapped or read.is_secondary:
#             continue
#         chrom = read.reference_name
#         if chrom not in ta_sites:
#             continue
#         # forward reads: insertion at read start
#         # reverse reads: insertion at read end
#         if read.is_reverse:
#             pos = read.reference_end  # 1-based end
#         else:
#             pos = read.reference_start + 1  # convert to 1-based
#         key = (chrom, pos)
#         counts[key] = counts.get(key, 0) + 1
#     bam.close()

#     # write WIG
#     wig_path = f"{WIGDIR}/{sample}.wig"
#     with open(wig_path, "w") as out:
#         current_chrom = None
#         for chrom, positions in ta_sites.items():
#             for pos in positions:
#                 if chrom != current_chrom:
#                     out.write(f"variableStep chrom={chrom}\n")
#                     current_chrom = chrom
#                 out.write(f"{pos}\t{counts.get((chrom, pos), 0)}\n")

#     total_insertions = sum(counts.values())
#     ta_with_insertions = len(counts)
#     print(f"  Done: {wig_path}")
#     print(f"  Total insertions: {total_insertions}")
#     print(f"  TA sites with insertions: {ta_with_insertions} / {total_sites}")

# print("All WIG files done")
# EOF

# # ── Step 4: Verify WIG files ──────────────────────────────────────
# echo "Verifying WIG files..."
# for SAMPLE in \
#     BHI_ERR1801012 BHI_ERR1801013 BHI_ERR1801014 \
#     HSerum_ERR1801009 HSerum_ERR1801010 HSerum_ERR1801011 \
#     Serum_ERR1801006 Serum_ERR1801007 Serum_ERR1801008; do
#         WIG=$WIGDIR/${SAMPLE}.wig
#         if [ ! -s $WIG ]; then
#             echo "ERROR: missing or empty WIG file: $WIG"
#             exit 1
#         fi
#         echo "  OK: $WIG ($(wc -l < $WIG) lines)"
# done
# echo "All WIG files verified"

mkdir -p $FILTERED

# filter WIG to active sites only
$PYTHON - <<EOF
WIGDIR = "$WIGDIR"
FILTERED = "$FILTERED"

samples = [
    "BHI_ERR1801012", "BHI_ERR1801013", "BHI_ERR1801014",
    "HSerum_ERR1801009", "HSerum_ERR1801010", "HSerum_ERR1801011",
    "Serum_ERR1801006", "Serum_ERR1801007", "Serum_ERR1801008"
]

# find all active sites across all samples
print("Finding active sites...")
active = set()
for sample in samples:
    current_chrom = None
    with open(f"{WIGDIR}/{sample}.wig") as f:
        for line in f:
            line = line.strip()
            if line.startswith("variableStep"):
                current_chrom = line.split("=")[1]
            else:
                pos, count = line.split("\t")
                if int(count) > 0:
                    active.add((current_chrom, int(pos)))
print(f"Active sites across all samples: {len(active)}")

# write filtered WIG files
for sample in samples:
    current_chrom = None
    last_written_chrom = None
    written = 0
    with open(f"{WIGDIR}/{sample}.wig") as f, \
         open(f"{FILTERED}/{sample}.wig", "w") as out:
        for line in f:
            line = line.strip()
            if line.startswith("variableStep"):
                current_chrom = line.split("=")[1]
            else:
                pos, count = line.split("\t")
                if (current_chrom, int(pos)) in active:
                    if current_chrom != last_written_chrom:
                        out.write(f"variableStep chrom={current_chrom}\n")
                        last_written_chrom = current_chrom
                    out.write(f"{pos}\t{count}\n")
                    written += 1
    print(f"  {sample}: {written} sites written")
print("Done")
EOF

# verify
echo "Lines in filtered WIG:"
wc -l $FILTERED/BHI_ERR1801012.wig


# run resampling with filtered WIGs
echo "Running resampling BHI vs HSerum..."
$PYTHON -m pytransit resampling \
    $FILTERED/BHI_ERR1801012.wig,$FILTERED/BHI_ERR1801013.wig,$FILTERED/BHI_ERR1801014.wig \
    $FILTERED/HSerum_ERR1801009.wig,$FILTERED/HSerum_ERR1801010.wig,$FILTERED/HSerum_ERR1801011.wig \
    $PROT_TABLE \
    $OUTDIR/BHI_vs_HSerum.txt \
    --s 10000 \
    || { echo "ERROR: resampling BHI vs HSerum failed"; exit 1; }
echo "BHI vs HSerum done"

echo "Running resampling BHI vs Serum..."
$PYTHON -m pytransit resampling \
    $FILTERED/BHI_ERR1801012.wig,$FILTERED/BHI_ERR1801013.wig,$FILTERED/BHI_ERR1801014.wig \
    $FILTERED/Serum_ERR1801006.wig,$FILTERED/Serum_ERR1801007.wig,$FILTERED/Serum_ERR1801008.wig \
    $PROT_TABLE \
    $OUTDIR/BHI_vs_Serum.txt \
    --s 10000 \
    || { echo "ERROR: resampling BHI vs Serum failed"; exit 1; }
echo "BHI vs Serum done"

echo "=== Results ==="
echo "Significant genes HSerum (q<0.05):"
awk 'NR>1 && $8 < 0.05' $OUTDIR/BHI_vs_HSerum.txt | wc -l
echo "Significant genes Serum (q<0.05):"
awk 'NR>1 && $8 < 0.05' $OUTDIR/BHI_vs_Serum.txt | wc -l
echo "Top 20 most depleted in Serum:"
awk 'NR>1' $OUTDIR/BHI_vs_Serum.txt | sort -k5 -n | head -20
