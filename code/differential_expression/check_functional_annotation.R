# ── Paths ─────────────────────────────────────────────────────────
OUTDIR_DE  <- "/home/mdha5802/genome_analysis/output/differential_expression/rnaseq"
OUTDIR_EGG <- "/proj/uppmax2026-1-61/nobackup/mdha5802/output/refined_annotation"

# ── Load eggNOG with correct headers ─────────────────────────────
# skip comment lines starting with ## but keep the header line starting with #query
egg_raw <- readLines(file.path(OUTDIR_EGG, "E745.emapper.annotations"))

# Find the header line (starts with #query)
header_line <- grep("^#query", egg_raw, value = TRUE)
cat("Header line found:\n", header_line, "\n\n")

# Read the actual data skipping all ## comment lines
egg <- read.table(
    file.path(OUTDIR_EGG, "E745.emapper.annotations"),
    sep = "\t",
    header = FALSE,
    comment.char = "",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    skip = sum(grepl("^##", egg_raw))   # skip ## comment lines only
)

# Set column names from header line
col_names <- unlist(strsplit(gsub("^#", "", header_line), "\t"))
colnames(egg) <- col_names

# Remove the header row if it got read as data
egg <- egg[egg$query != "query", ]

cat("Dimensions:", nrow(egg), "genes,", ncol(egg), "columns\n")
cat("Column names:\n")
print(colnames(egg))
cat("\nFirst 3 rows:\n")
print(head(egg[, c("query", "Preferred_name", "COG_category", 
                    "Description", "KEGG_Pathway")], 3))