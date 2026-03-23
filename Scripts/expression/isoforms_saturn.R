#!/usr/bin/env Rscript
# Differential transcript usage between standard and depleted
# Every RNA library contains a matched standard and depleted hence the pair term

suppressPackageStartupMessages({
  library(tximport)
  library(edgeR)
  library(SummarizedExperiment)
  library(readr)
  library(dplyr)
  library(stringr)
  library(limma)
  library(BiocParallel)
  library(satuRn, lib.loc="~/Rpackages")
})

# -------------------------
# 1) Files & sample table
# -------------------------
files <- list.files(pattern = "\\.isoforms\\.results$", full.names = TRUE)
stopifnot(length(files) > 0)

samples <- tibble(
  file       = files,
  sample_id  = basename(files) |> str_remove("\\.isoforms\\.results$"),
  condition  = if_else(str_detect(sample_id, "-WM$"), "WM", "D"),
  pair       = str_replace(sample_id, "-(WM|D)$", "")  
) |>
  mutate(
    condition = factor(condition, levels = c("WM","D")),  
    pair      = factor(make.names(as.character(pair)))
  )

# Sanity: every pair has both WM and D
tbl_pair <- with(samples, table(pair, condition))
stopifnot(all(tbl_pair[, "WM"] >= 1 & tbl_pair[, "D"] >= 1))

# -------------------------
# 2) Import transcript counts
# -------------------------
txi <- tximport(files = samples$file,
                type  = "rsem",
                txIn  = TRUE,
                txOut = TRUE,
                countsFromAbundance = "no")   

cts <- round(txi$counts)
stopifnot(ncol(cts) == nrow(samples))
colnames(cts) <- samples$sample_id

# -------------------------
# 3) Build transcript -> gene map
# -------------------------
message("reading in files with read_tsv")
tx2gene <- bind_rows(lapply(samples$file, function(f)
  read_tsv(f, col_types = cols(
    transcript_id = col_character(),
    gene_id       = col_character(),
    .default      = col_skip()
  )) |> select(isoform_id = transcript_id, gene_id)
)) |> distinct()

# Keep only transcripts present in cts and ensure matching order
tx2gene <- tx2gene[match(rownames(cts), tx2gene$isoform_id), ]
stopifnot(all(tx2gene$isoform_id == rownames(cts)))

# -------------------------
# 4) SummarizedExperiment
# -------------------------
se <- SummarizedExperiment(
  assays  = list(counts = cts),
  colData = DataFrame(samples),
  rowData = DataFrame(tx2gene)    
)
rownames(se) <- rownames(cts)

# -------------------------
# 5) Filtering
# -------------------------
keep_tx <- filterByExpr(assay(se, "counts"), group = se$condition)
se <- se[keep_tx, ]

# Enforce ≥2 transcripts per gene
tab <- table(rowData(se)$gene_id)
se  <- se[rowData(se)$gene_id %in% names(tab[tab >= 2]), ]

message("Dims after filtering: ", paste(dim(se), collapse = " x "))

# -------------------------
# 6) Fit satuRn (paired design)
# -------------------------
bpp <- MulticoreParam(workers = max(1, parallel::detectCores() - 1))

se <- satuRn::fitDTU(
  object   = se,
  formula  = ~ pair + condition,
  parallel = TRUE,
  BPPARAM  = bpp,
  verbose  = TRUE
)

# -------------------------
# 7) Build contrast & CHUNKED testing to avoid locfdr path
# -------------------------
design <- model.matrix(~ pair + condition, data = as.data.frame(colData(se)))
coef_names <- colnames(design)
stopifnot("conditionD" %in% coef_names)

L <- matrix(0, nrow = length(coef_names), ncol = 1,
            dimnames = list(coef_names, "D_vs_WM"))
L["conditionD", 1] <- 1

# Split into blocks smaller than ~500 to bypass empirical-null step
chunk_size <- 400L
idx_list <- split(seq_len(nrow(se)), ceiling(seq_len(nrow(se))/chunk_size))

get_block <- function(i) {
  sei <- se[idx_list[[i]], ]
  sei <- satuRn::testDTU(
    object    = sei,
    contrasts = L,
    diagplot1 = FALSE,
    diagplot2 = FALSE,
    sort      = FALSE
  )
  out <- as.data.frame(rowData(sei)[["fitDTUResult_D_vs_WM"]])
  out$isoform_id <- rownames(sei)
  # carry gene_id directly from this block's rowData
  out$gene_id <- rowData(sei)$gene_id
  out
}

message("Testing in ", length(idx_list), " chunks of ≤ ", chunk_size, " transcripts...")
res_list <- lapply(seq_along(idx_list), get_block)
res_tx   <- bind_rows(res_list)

# -------------------------
# 8) Adjust p-values across ALL transcripts
# -------------------------
# Ensure numeric p-values
res_tx$pval <- as.numeric(res_tx$pval)
res_tx$padj <- p.adjust(res_tx$pval, method = "BH")

# Reorder columns nicely
res_tx <- res_tx |>
  relocate(isoform_id, gene_id) |>
  arrange(padj)

# Save transcript-level results
write.csv(res_tx, "satuRn_DTU_transcripts_D_vs_WM.csv", row.names = FALSE)

# -------------------------
# 9) Add Δusage (D - WM) effect size (descriptive)
# -------------------------
cts_cur <- assay(se, "counts")
gene_id_vec <- rowData(se)$gene_id

gene_totals <- rowsum(cts_cur, group = gene_id_vec)
denom <- gene_totals[ match(gene_id_vec, rownames(gene_totals)), , drop = FALSE ]
usage <- cts_cur / pmax(denom, 1)  

wm_idx <- se$condition == "WM"
d_idx  <- se$condition == "D"

mean_WM <- rowMeans(usage[, wm_idx, drop = FALSE])
mean_D  <- rowMeans(usage[, d_idx,  drop = FALSE])

delta_df <- tibble(
  isoform_id     = rownames(se),
  mean_usage_WM  = as.numeric(mean_WM),
  mean_usage_D   = as.numeric(mean_D),
  delta_usage    = mean_usage_D - mean_usage_WM
)

res_tx_du <- res_tx |>
  left_join(delta_df, by = "isoform_id") |>
  arrange(padj)

write.csv(res_tx_du, "satuRn_DTU_transcripts_D_vs_WM_with_deltaUsage.csv", row.names = FALSE)

# -------------------------
# 10) (Optional) Per-gene Simes aggregation (quick gene-level DTU)
# -------------------------
simes <- function(p) { p <- sort(p); m <- length(p); min(1, min(m * p / seq_along(p))) }
res_gene <- res_tx |>
  group_by(gene_id) |>
  summarize(p_simes = simes(pval), .groups = "drop") |>
  mutate(padj_simes = p.adjust(p_simes, "BH")) |>
  arrange(padj_simes)

write.csv(res_gene, "satuRn_DTU_genes_Simes.csv", row.names = FALSE)

message("Done. Wrote:",
        "\n  - satuRn_DTU_transcripts_D_vs_WM.csv",
        "\n  - satuRn_DTU_transcripts_D_vs_WM_with_deltaUsage.csv",
        "\n  - satuRn_DTU_genes_Simes.csv")


