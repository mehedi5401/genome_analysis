# 1. Feature type breakdown
grep -v "^#" E745.gff | awk '{print $3}' | sort | uniq -c | sort -rn

# 2. Hypothetical protein count and percentage
total=$(grep -c "^>" E745.faa)
hypo=$(grep -c "hypothetical protein" E745.faa)
echo "Total CDSs: $total"
echo "Hypothetical proteins: $hypo"
awk "BEGIN {printf \"Hypothetical %%: %.1f\n\", ($hypo/$total)*100}"

# 3. Locus tag format
grep -o 'locus_tag=[^;]*' E745.gff | head -5

# 4. Prokka version and command used (if log exists)
cat E745.log 2>/dev/null | head -20 || echo "No log found"

# 5. Contig count fed to Prokka
grep -c "^>" E745.contigs.fasta