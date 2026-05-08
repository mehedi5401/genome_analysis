.libPaths('/home/mdha5802/R/library')

# ── Set CRAN mirror ──────────────────────────────────────────────
options(repos = c(CRAN = 'https://cloud.r-project.org'))

# ── Install packages if not present ─────────────────────────────
if (!require('BiocManager', quietly = FALSE))
    install.packages('BiocManager', lib = '/home/mdha5802/R/library')

if (!require('DESeq2', quietly = FALSE))
    BiocManager::install('DESeq2', lib = '/home/mdha5802/R/library')

if (!require('ggplot2', quietly = FALSE))
    install.packages('ggplot2', lib = '/home/mdha5802/R/library')

library(DESeq2)
library(ggplot2)

# ── Paths ────────────────────────────────────────────────────────
COUNTDIR <- "/home/mdha5802/genome_analysis/output/counting/rnaseq"
OUTDIR   <- "/home/mdha5802/genome_analysis/output/differential_expression/rnaseq"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ── Load count files ─────────────────────────────────────────────
samples <- data.frame(
    file      = c("ERR1797972_counts.txt",
                  "ERR1797973_counts.txt",
                  "ERR1797974_counts.txt",
                  "ERR1797969_counts.txt",
                  "ERR1797970_counts.txt",
                  "ERR1797971_counts.txt"),
    condition = factor(c("BHI", "BHI", "BHI",
                         "Serum", "Serum", "Serum")),
    row.names = c("BHI_1", "BHI_2", "BHI_3",
                  "Serum_1", "Serum_2", "Serum_3")
)

# Read all count files into one matrix
counts <- do.call(cbind, lapply(samples$file, function(f) {
    df <- read.table(file.path(COUNTDIR, f),
                     header = FALSE,
                     row.names = 1,
                     sep = "\t")
    df[, 1]
}))
rownames(counts) <- rownames(read.table(
    file.path(COUNTDIR, samples$file[1]),
    header = FALSE, row.names = 1, sep = "\t"))
colnames(counts) <- rownames(samples)

# Remove HTSeq summary lines (start with __)
counts <- counts[!grepl("^__", rownames(counts)), ]

cat("Genes loaded:", nrow(counts), "\n")
cat("Samples loaded:", ncol(counts), "\n")

# ── Create DESeq2 object ─────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = samples,
    design    = ~ condition
)

# Set BHI as reference condition
dds$condition <- relevel(dds$condition, ref = "BHI")

# Filter low count genes — keep genes with at least 10 reads total
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes after filtering:", nrow(dds), "\n")

# ── Run DESeq2 ───────────────────────────────────────────────────
cat("Running DESeq2...\n")
dds <- DESeq(dds)

# ── Extract results ──────────────────────────────────────────────
# Positive log2FC = higher in Serum
# Negative log2FC = lower in Serum
res <- results(dds, contrast = c("condition", "Serum", "BHI"))
res <- res[order(res$padj), ]

cat("Total DE genes (padj < 0.05):",
    sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in Serum (padj < 0.05, log2FC > 1):",
    sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("Downregulated in Serum (padj < 0.05, log2FC < -1):",
    sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE), "\n")

# Save results
write.csv(as.data.frame(res),
          file.path(OUTDIR, "serum_vs_BHI_results.csv"))
cat("Results saved to:", file.path(OUTDIR, "serum_vs_BHI_results.csv"), "\n")

# ── Normalised counts ────────────────────────────────────────────
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts,
          file.path(OUTDIR, "normalised_counts.csv"))

# ── PCA plot ─────────────────────────────────────────────────────
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.8, size = 3) +
    xlab(paste0("PC1: ", pca_var[1], "% variance")) +
    ylab(paste0("PC2: ", pca_var[2], "% variance")) +
    theme_bw()
ggsave(file.path(OUTDIR, "pca_plot.pdf"), width = 6, height = 5)
cat("PCA plot saved\n")

# ── Volcano plot ─────────────────────────────────────────────────
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$sig  <- "not significant"
res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "up in serum"
res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "down in serum"

ggplot(res_df, aes(x = log2FoldChange,
                   y = -log10(padj),
                   color = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c(
        "not significant" = "grey60",
        "up in serum"     = "#D85A30",
        "down in serum"   = "#185FA5")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = "Serum vs BHI",
         x = "Log2 fold change",
         y = "-Log10 adjusted p-value",
         color = "") +
    theme_bw()
ggsave(file.path(OUTDIR, "volcano_plot.pdf"), width = 8, height = 6)
cat("Volcano plot saved\n")

cat("DESeq2 analysis complete\n")