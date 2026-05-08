.libPaths('/home/mdha5802/R/library')
# library(DESeq2)

# ── Load results ─────────────────────────────────────────────────
OUTDIR <- "/home/mdha5802/genome_analysis/output/differential_expression/rnaseq"
res_df <- read.csv(file.path(OUTDIR, "serum_vs_BHI_results.csv"), row.names = 1)

# ── Parse GFF annotation ─────────────────────────────────────────
gff_file <- "/home/mdha5802/genome_analysis/output/annotation/E745_clean.gff"
gff <- read.table(gff_file, sep = "\t", quote = "", comment.char = "#",
                  col.names = c("seqid","source","type","start","end",
                                "score","strand","phase","attributes"))

# Keep only gene/CDS features
gff <- gff[gff$type %in% c("gene", "CDS"), ]

# ── Extract attributes using gsub (more robust than regmatches) ───
extract_attr <- function(attributes, key) {
    pattern <- paste0(".*", key, "=([^;]+).*")
    result <- ifelse(grepl(paste0(key, "="), attributes),
                     gsub(pattern, "\\1", attributes),
                     NA_character_)
    return(as.character(result))
}

gff$locus_tag <- extract_attr(gff$attributes, "locus_tag")
gff$gene_name <- extract_attr(gff$attributes, "gene")
gff$product   <- extract_attr(gff$attributes, "product")

# Keep one row per locus_tag
gff_annot <- gff[!is.na(gff$locus_tag), c("locus_tag","gene_name","product","strand")]
gff_annot <- gff_annot[!duplicated(gff_annot$locus_tag), ]

# Ensure all columns are plain character vectors
gff_annot$locus_tag <- as.character(gff_annot$locus_tag)
gff_annot$gene_name <- as.character(gff_annot$gene_name)
gff_annot$product   <- as.character(gff_annot$product)
gff_annot$strand    <- as.character(gff_annot$strand)

# ── Merge annotation into results ────────────────────────────────
res_df$locus_tag <- rownames(res_df)
res_annotated <- merge(res_df, gff_annot, by = "locus_tag", all.x = TRUE)
rownames(res_annotated) <- res_annotated$locus_tag

# ── Verify no list columns remain ────────────────────────────────
list_cols <- sapply(res_annotated, is.list)
if (any(list_cols)) {
    cat("Converting list columns:", names(which(list_cols)), "\n")
    res_annotated[list_cols] <- lapply(res_annotated[list_cols], 
                                        function(x) sapply(x, paste, collapse = ";"))
}

# ── Save full annotated results ───────────────────────────────────
write.csv(res_annotated,
          file.path(OUTDIR, "serum_vs_BHI_results_annotated.csv"),
          row.names = FALSE)
cat("Annotated results saved\n")

# ── Top 10 most significant with gene names ───────────────────────
cat("\nTop 10 most significant DE genes:\n")
top10_sig <- head(res_annotated[order(res_annotated$padj),
                   c("locus_tag","gene_name","product","baseMean",
                     "log2FoldChange","padj")], 10)
print(top10_sig)

# ── Top 10 most upregulated with gene names ───────────────────────
cat("\nTop 10 most upregulated in Serum:\n")
top10_up <- head(res_annotated[order(-res_annotated$log2FoldChange),
                  c("locus_tag","gene_name","product","baseMean",
                    "log2FoldChange","padj")], 10)
print(top10_up)

# ── Summary counts ────────────────────────────────────────────────
cat("\nTotal DE genes (padj < 0.05):",
    sum(res_annotated$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in Serum (padj < 0.05, log2FC > 1):",
    sum(res_annotated$padj < 0.05 & res_annotated$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("Downregulated in Serum (padj < 0.05, log2FC < -1):",
    sum(res_annotated$padj < 0.05 & res_annotated$log2FoldChange < -1, na.rm = TRUE), "\n")