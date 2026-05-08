# Find locus tags for key genes in Prokka GFF
grep -i "pyrD\|purD\|manY\|pyrK\|pyrF\|purH" \
    /home/mdha5802/genome_analysis/output/annotation/E745_clean.gff | \
    grep "CDS" | \
    awk '{print $ $9}' | \
    tr ';' '\n' | \
    grep "ID=\|product="