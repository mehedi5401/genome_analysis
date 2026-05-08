.libPaths('/home/mdha5802/R/library')

# ── Set CRAN mirror ──────────────────────────────────────────────
options(repos = c(CRAN = 'https://cloud.r-project.org'))

# ── Install packages if not present ─────────────────────────────
if (!require('ggplot2', quietly = FALSE))
    install.packages('ggplot2', lib = '/home/mdha5802/R/library')

library(ggplot2)


# ── Paths ─────────────────────────────────────────────────────────
OUTDIR_DE  <- "/home/mdha5802/genome_analysis/output/differential_expression/rnaseq"
OUTDIR_EGG <- "/proj/uppmax2026-1-61/nobackup/mdha5802/output/refined_annotation"
OUTDIR_OUT <- "/home/mdha5802/genome_analysis/output/differential_expression/rnaseq"

# ── Load eggNOG ───────────────────────────────────────────────────
egg_raw <- readLines(file.path(OUTDIR_EGG, "E745.emapper.annotations"))
header_line <- grep("^#query", egg_raw, value = TRUE)

egg <- read.table(
    file.path(OUTDIR_EGG, "E745.emapper.annotations"),
    sep = "\t", header = FALSE, comment.char = "",
    quote = "", fill = TRUE, stringsAsFactors = FALSE,
    skip = sum(grepl("^##", egg_raw)))

col_names <- unlist(strsplit(gsub("^#", "", header_line), "\t"))
colnames(egg) <- col_names
egg <- egg[egg$query != "query", ]
egg[egg == "-"] <- NA

# ── Load DESeq2 results ───────────────────────────────────────────
res <- read.csv(file.path(OUTDIR_DE, "serum_vs_BHI_results_annotated.csv"))
cat("DESeq2 results loaded:", nrow(res), "genes\n")

# ── Merge ─────────────────────────────────────────────────────────
res_full <- merge(res,
                  egg[, c("query", "Preferred_name", "COG_category",
                           "Description", "KEGG_Pathway", "KEGG_Module",
                           "KEGG_ko", "GOs", "PFAMs")],
                  by.x = "locus_tag", by.y = "query",
                  all.x = TRUE)

cat("Merged dimensions:", nrow(res_full), "genes\n")
cat("Genes with COG annotation:", sum(!is.na(res_full$COG_category)), "\n")
cat("Genes with KEGG pathway:",   sum(!is.na(res_full$KEGG_Pathway)), "\n")

# ── Save full merged results ──────────────────────────────────────
write.csv(res_full,
          file.path(OUTDIR_OUT, "serum_vs_BHI_full_annotated.csv"),
          row.names = FALSE)
cat("Saved: serum_vs_BHI_full_annotated.csv\n")

# ── DE gene subsets ───────────────────────────────────────────────
de_all  <- res_full[!is.na(res_full$padj) & res_full$padj < 0.05, ]
de_up   <- de_all[de_all$log2FoldChange >  1, ]
de_down <- de_all[de_all$log2FoldChange < -1, ]

cat("\nTotal DE genes:         ", nrow(de_all),  "\n")
cat("Upregulated in Serum:   ", nrow(de_up),   "\n")
cat("Downregulated in Serum: ", nrow(de_down), "\n")

# ── COG category reference ────────────────────────────────────────
cog_ref <- c(
    J = "Translation",
    K = "Transcription",
    L = "Replication & repair",
    M = "Cell wall/membrane biogenesis",
    N = "Cell motility",
    O = "Post-translational modification",
    P = "Inorganic ion transport & metabolism",
    Q = "Secondary metabolite biosynthesis",
    T = "Signal transduction",
    U = "Intracellular trafficking & secretion",
    V = "Defense mechanisms",
    W = "Extracellular structures",
    C = "Energy production & conversion",
    D = "Cell cycle control & division",
    E = "Amino acid transport & metabolism",
    F = "Nucleotide transport & metabolism",
    G = "Carbohydrate transport & metabolism",
    H = "Coenzyme transport & metabolism",
    I = "Lipid transport & metabolism",
    R = "General function prediction only",
    S = "Function unknown"
)

# ── COG summary function ──────────────────────────────────────────
summarise_cog <- function(df, label) {
    cogs <- unlist(strsplit(na.omit(df$COG_category), ""))
    cogs <- cogs[cogs %in% names(cog_ref)]
    tbl  <- sort(table(cogs), decreasing = TRUE)
    tbl_df <- data.frame(
        COG      = names(tbl),
        Count    = as.integer(tbl),
        Function = cog_ref[names(tbl)]
    )
    cat("\n── COG categories —", label, "──\n")
    print(tbl_df)
    return(tbl_df)
}

cog_all  <- summarise_cog(de_all,  "All DE genes")
cog_up   <- summarise_cog(de_up,   "Upregulated in Serum")
cog_down <- summarise_cog(de_down, "Downregulated in Serum")

# ── COG bar plot ──────────────────────────────────────────────────
library(ggplot2)

cog_up$direction   <- "Up in Serum"
cog_down$direction <- "Down in Serum"
cog_plot_df <- rbind(cog_up, cog_down)

ggplot(cog_plot_df, aes(x = reorder(COG, Count),
                         y = Count,
                         fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Up in Serum"   = "#D85A30",
                                  "Down in Serum" = "#185FA5")) +
    coord_flip() +
    labs(title = "COG Categories — DE Genes Serum vs BHI",
         x = "COG Category",
         y = "Number of genes",
         fill = "") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))

ggsave(file.path(OUTDIR_OUT, "COG_barplot.pdf"), width = 8, height = 6)
cat("COG barplot saved\n")

# ── Top DE genes with full annotation ────────────────────────────
cat("\nTop 10 upregulated with functional annotation:\n")
print(head(de_up[order(-de_up$log2FoldChange),
                  c("locus_tag", "Preferred_name", "Description",
                    "COG_category", "log2FoldChange", "padj")], 10))

cat("\nTop 10 downregulated with functional annotation:\n")
print(head(de_down[order(de_down$log2FoldChange),
                    c("locus_tag", "Preferred_name", "Description",
                      "COG_category", "log2FoldChange", "padj")], 10))