# Remove FASTA sequence from GFF
awk '/^##FASTA/{exit} {print}' \
    /home/mdha5802/genome_analysis/output/annotation/E745.gff \
    > /home/mdha5802/genome_analysis/output/annotation/E745_clean.gff

# Verify — should show only GFF lines, no FASTA
tail -5 /home/mdha5802/genome_analysis/output/annotation/E745_clean.gff
wc -l /home/mdha5802/genome_analysis/output/annotation/E745_clean.gff