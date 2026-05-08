.libPaths('/home/mdha5802/R/library')
# ── Load annotated results ────────────────────────────────────────
OUTDIR <- "/home/mdha5802/genome_analysis/output/differential_expression/rnaseq"
res_annotated <- read.csv(file.path(OUTDIR, "serum_vs_BHI_results_annotated.csv"))

# ── Target genes ──────────────────────────────────────────────────
target_genes <- c("pyrD", "purD", "manY", "pyrK", "pyrF", "purH", "pyrK_2", "manY_2")

# ── Search by gene_name column ────────────────────────────────────
hits <- res_annotated[
    !is.na(res_annotated$gene_name) &
    tolower(res_annotated$gene_name) %in% tolower(target_genes), ]

# ── Also search by product column in case gene_name is missing ────
hits_product <- res_annotated[
    !is.na(res_annotated$product) &
    grepl(paste(target_genes, collapse = "|"), 
          res_annotated$product, ignore.case = TRUE), ]

# ── Combine both hits ─────────────────────────────────────────────
hits_combined <- unique(rbind(hits, hits_product))

# ── Print results ─────────────────────────────────────────────────
cat("Found", nrow(hits_combined), "matching genes\n\n")

print(hits_combined[, c("locus_tag", "gene_name", "product",
                         "baseMean", "log2FoldChange", 
                         "lfcSE", "pvalue", "padj")])

# ── Flag significance ─────────────────────────────────────────────
hits_combined$status <- "not significant"
hits_combined$status[hits_combined$padj < 0.05 & 
                     hits_combined$log2FoldChange > 1]  <- "up in serum"
hits_combined$status[hits_combined$padj < 0.05 & 
                     hits_combined$log2FoldChange < -1] <- "down in serum"

cat("\nExpression status:\n")
print(hits_combined[, c("locus_tag", "gene_name", "log2FoldChange", 
                         "padj", "status")])

# ── Save to file ──────────────────────────────────────────────────
write.csv(hits_combined, 
          file.path(OUTDIR, "target_genes_expression.csv"),
          row.names = FALSE)
cat("\nSaved to target_genes_expression.csv\n")