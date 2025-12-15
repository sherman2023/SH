#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# Allow running from anywhere by taking paths from environment variables.
samples_tsv <- Sys.getenv("SAMPLES_TSV", "config/samples.tsv")
ref_dir     <- Sys.getenv("REF_DIR", "ref")
results_dir <- Sys.getenv("RESULTS_DIR", "results")

dir.create(file.path(results_dir, "deseq2"), recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Sample metadata
# ----------------------------
samples <- read_tsv(samples_tsv, show_col_types = FALSE) %>%
  mutate(cell = factor(cell),
         dex  = factor(dex, levels = c("untrt","trt")))

# ----------------------------
# Approach A: featureCounts -> DESeq2
# ----------------------------
fc_path <- file.path(results_dir, "approachA_featurecounts", "featureCounts.counts.tsv")
stopifnot(file.exists(fc_path))

fc <- read_tsv(fc_path, comment = "#", show_col_types = FALSE)

# featureCounts format: annotation columns then count columns begin at column 7
count_cols <- colnames(fc)[7:ncol(fc)]
count_mat <- as.matrix(fc[, count_cols])
rownames(count_mat) <- fc$Geneid

# Clean sample names: remove paths and STAR BAM suffix
colnames(count_mat) <- gsub(".*/", "", colnames(count_mat))
colnames(count_mat) <- gsub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(count_mat))

# Reorder columns to match samples.tsv
count_mat <- count_mat[, samples$run]

ddsA <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData   = as.data.frame(samples),
                              design    = ~ cell + dex)
ddsA <- DESeq(ddsA)
resA <- results(ddsA, contrast = c("dex","trt","untrt")) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

write_tsv(resA, file.path(results_dir, "deseq2", "resA_featureCounts_DESeq2.tsv"))

# ----------------------------
# Approach B: Salmon -> tximport -> DESeq2
# ----------------------------
tx2gene_path <- file.path(ref_dir, "tx2gene.tsv")
stopifnot(file.exists(tx2gene_path))
tx2gene <- read_tsv(tx2gene_path, col_names = c("tx","gene"), show_col_types = FALSE)

files <- file.path(results_dir, "approachB_salmon", samples$run, "quant.sf")
names(files) <- samples$run
stopifnot(all(file.exists(files)))

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

ddsB <- DESeqDataSetFromTximport(txi,
                                colData = as.data.frame(samples),
                                design  = ~ cell + dex)
ddsB <- DESeq(ddsB)
resB <- results(ddsB, contrast = c("dex","trt","untrt")) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

write_tsv(resB, file.path(results_dir, "deseq2", "resB_salmon_tximport_DESeq2.tsv"))

# ----------------------------
# Compare results (A vs B)
# ----------------------------
comp <- resA %>%
  select(gene_id,
         log2FC_A = log2FoldChange,
         padj_A   = padj) %>%
  inner_join(resB %>% select(gene_id,
                             log2FC_B = log2FoldChange,
                             padj_B   = padj),
             by = "gene_id")

sigA <- comp %>% filter(!is.na(padj_A), padj_A < 0.05) %>% pull(gene_id) %>% unique()
sigB <- comp %>% filter(!is.na(padj_B), padj_B < 0.05) %>% pull(gene_id) %>% unique()

overlap <- length(intersect(sigA, sigB))
jaccard <- ifelse(length(union(sigA, sigB)) == 0, NA, overlap / length(union(sigA, sigB)))

stats <- tibble(
  n_genes_compared = nrow(comp),
  n_sig_A_padj_0.05 = length(sigA),
  n_sig_B_padj_0.05 = length(sigB),
  n_overlap = overlap,
  jaccard = jaccard,
  cor_log2FC = cor(comp$log2FC_A, comp$log2FC_B, use = "pairwise.complete.obs")
)

write_tsv(stats, file.path(results_dir, "deseq2", "AvsB_summary_stats.tsv"))

p_scatter <- ggplot(comp, aes(x = log2FC_A, y = log2FC_B)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  labs(title = "A vs B: log2 fold-changes (DESeq2)",
       x = "log2FC (STAR+featureCounts)",
       y = "log2FC (Salmon+tximport)") +
  theme_minimal(base_size = 12)

ggsave(file.path(results_dir, "deseq2", "AvsB_log2FC_scatter.png"),
       p_scatter, width = 6, height = 5, dpi = 300)

volcano <- function(df, title, out_png) {
  df2 <- as_tibble(df) %>%
    mutate(minusLog10Padj = -log10(padj)) %>%
    mutate(minusLog10Padj = ifelse(is.infinite(minusLog10Padj), NA, minusLog10Padj))
  p <- ggplot(df2, aes(x = log2FoldChange, y = minusLog10Padj)) +
    geom_point(alpha = 0.25, size = 0.8) +
    labs(title = title, x = "log2 fold-change", y = "-log10(adjusted p-value)") +
    theme_minimal(base_size = 12)
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
}

volcano(as.data.frame(results(ddsA, contrast=c("dex","trt","untrt"))),
        "Volcano: STAR+featureCounts → DESeq2",
        file.path(results_dir, "deseq2", "volcano_A.png"))

volcano(as.data.frame(results(ddsB, contrast=c("dex","trt","untrt"))),
        "Volcano: Salmon+tximport → DESeq2",
        file.path(results_dir, "deseq2", "volcano_B.png"))

sink(file.path(results_dir, "deseq2", "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. See: ", file.path(results_dir, "deseq2"))
