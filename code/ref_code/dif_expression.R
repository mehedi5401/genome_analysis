# Install DESeq2 if needed (only once)
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)

# ── 1. Load count data ──────────────────────────────────────────────

counts_raw <- read.table(
    "1_Zhang_2017/rnaseq/counts/all_samples_counts.txt",
    header = TRUE,
    skip = 1,      # featureCounts adds a header line to skip
    row.names = 1
)

# Keep only the count columns (remove chromosome, start, end, strand, length)
counts <- counts_raw[, 6:ncol(counts_raw)]

# Clean up column names to just sample names
colnames(counts) <- c("BHI_1", "BHI_2", "BHI_3",
                       "Serum_1", "Serum_2", "Serum_3")

# ── 2. Define sample information ────────────────────────────────────

sample_info <- data.frame(
    condition = factor(c("BHI", "BHI", "BHI",
                         "Serum", "Serum", "Serum")),
    row.names = colnames(counts)
)

# ── 3. Create DESeq2 object ─────────────────────────────────────────

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sample_info,
    design    = ~ condition
)

# ── 4. Run DESeq2 ───────────────────────────────────────────────────

dds <- DESeq(dds)

# ── 5. Extract results ──────────────────────────────────────────────
# Serum vs BHI: positive = higher in serum, negative = lower in serum

results <- results(dds,
                   contrast = c("condition", "Serum", "BHI"))

results_df <- as.data.frame(results)
results_df <- results_df[order(results_df$padj), ]  # sort by significance

# Save results
dir.create("1_Zhang_2017/rnaseq/deseq2", showWarnings = FALSE)

write.csv(results_df,
          "1_Zhang_2017/rnaseq/deseq2/serum_vs_BHI_results.csv")

# ── 6. Summary ──────────────────────────────────────────────────────

cat("Upregulated in serum (padj < 0.05, log2FC > 1):",
    sum(results_df$padj < 0.05 & results_df$log2FoldChange > 1,
        na.rm = TRUE), "\n")

cat("Downregulated in serum (padj < 0.05, log2FC < -1):",
    sum(results_df$padj < 0.05 & results_df$log2FoldChange < -1,
        na.rm = TRUE), "\n")

# ── 7. Volcano plot ─────────────────────────────────────────────────

results_df$significant <- "Not significant"
results_df$significant[results_df$padj < 0.05 &
                        results_df$log2FoldChange > 1] <- "Upregulated in serum"
results_df$significant[results_df$padj < 0.05 &
                        results_df$log2FoldChange < -1] <- "Downregulated in serum"

ggplot(results_df, aes(x = log2FoldChange,
                        y = -log10(padj),
                        color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("grey", "red", "blue")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = "Serum vs BHI — Differential Expression",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_bw()

ggsave("1_Zhang_2017/rnaseq/deseq2/volcano_plot.pdf",
       width = 8, height = 6)

cat("DESeq2 analysis complete\n")