#!/usr/bin/env Rscript
# DESeq2 over isoforms between standard and depleted 
# Every RNA library contains a matched standard and depleted hence the pair term

library(tidyverse)
library(tximport)
library(DESeq2)

# 1) Point to your RSEM isoform files
files <- list.files(pattern = "\\.isoforms\\.results$", full.names = TRUE)
stopifnot(length(files) > 0)

# 2) Build sample table from filenames like "1-3-WM.isoforms.results"
samples <- tibble(
  file      = files,
  sample    = basename(files) |> str_remove("\\.isoforms\\.results$"),
  condition = if_else(str_detect(sample, "-WM$"), "WM", "D"),
  pair      = str_replace(sample, "-(WM|D)$", "")   
)

# Build samples (ensure factors)
samples <- samples %>%
  mutate(condition = factor(condition, levels = c("WM","D")),
         pair      = factor(pair))

txi <- tximport(files = samples$file,
                type  = "rsem",
                txIn  = TRUE,         
                txOut = TRUE,         
                countsFromAbundance = "lengthScaledTPM")

# sanitize lengths (RSEM can give 0 effective lengths)
txi$length[!is.finite(txi$length) | txi$length <= 0] <- 1

#convert '-' to safe '_'
samples$pair <- gsub("[^A-Za-z0-9_.]", "_", samples$pair)

dds <- DESeqDataSetFromTximport(txi,
                                colData = as.data.frame(samples),
                                design = ~ pair + condition)

# need 1 counts (TPM1) in atleast 10 samples 
keep <- rowSums(counts(dds) >= 1) >= ceiling(0.1 * ncol(dds))
dds  <- dds[keep,]

# 5) Fit model and test D vs WM
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition","D","WM"))

# --- add group means (normalized + model-based, pair-adjusted) ---

# 1) Empirical group means from normalized counts
norm <- counts(dds, normalized = TRUE)
wm_idx <- colData(dds)$condition == "WM"
d_idx  <- colData(dds)$condition == "D"

mean_WM_norm <- rowMeans(norm[, wm_idx, drop = FALSE])
mean_D_norm  <- rowMeans(norm[, d_idx,  drop = FALSE])

# 2) Model-based, pair-adjusted means from fitted values (NB-GLM with ~ pair + condition)
#    - average fitted means within each pair/group, then average across pairs.
stopifnot("mu" %in% names(assays(dds)))  
mu   <- assays(dds)[["mu"]]              
pair <- as.character(colData(dds)$pair)

pairs_complete <- intersect(unique(pair[wm_idx]), unique(pair[d_idx]))
stopifnot(length(pairs_complete) > 0)

get_means_by_pair <- function(mu, group_idx, pair, pairs_complete) {
  sapply(pairs_complete, function(p) {
    idx <- which(group_idx & pair == p)
    if (length(idx) == 0) rep(NA_real_, nrow(mu)) else rowMeans(mu[, idx, drop = FALSE])
  })
}

mu_WM_by_pair <- get_means_by_pair(mu, wm_idx, pair, pairs_complete)
mu_D_by_pair  <- get_means_by_pair(mu, d_idx,  pair, pairs_complete)

mean_WM_model <- rowMeans(mu_WM_by_pair, na.rm = TRUE)
mean_D_model  <- rowMeans(mu_D_by_pair,  na.rm = TRUE)

# 3) Bind onto DESeq2 results and save
res_tbl <- as_tibble(res, rownames = "transcript_id") %>%
  mutate(
    mean_WM_norm  = mean_WM_norm[transcript_id],
    mean_D_norm   = mean_D_norm[transcript_id],
    mean_WM_model = mean_WM_model[transcript_id],
    mean_D_model  = mean_D_model[transcript_id]
  ) %>%
  arrange(padj)

write.csv(res_tbl, "DESeq2_D_vs_WM_paired_results_isoforms_with_means.csv", row.names = FALSE)

