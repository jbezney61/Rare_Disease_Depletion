#!/usr/bin/env Rscript
# used to annoate leafcuter junctions

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---- CLI ----
option_list <- list(
  make_option(c("-a","--annotation_code"), type="character",
              help="Prefix path for Leafcutter annotation set, e.g. annotation_codes/gencode_hg38/gencode_hg38"),
  make_option(c("-o","--output"), type="character", default=NULL,
              help="Output text file (default: <input>.annotated.txt)")
)
opt <- parse_args(OptionParser(
  usage="%prog -a <annotation_code_prefix> -o <out.txt> <perind_numers.counts.gz>",
  option_list=option_list
), positional_arguments = 1)

counts_gz <- opt$args[[1]]
if (is.null(opt$options$annotation_code))
  stop("Please provide -a/--annotation_code prefix (e.g. annotation_codes/gencode_hg38/gencode_hg38)")
anno_prefix <- opt$options$annotation_code

# Output name: plain text (no gzip)
if (is.null(opt$options$output)) {
  base  <- sub("\\.gz$", "", counts_gz, perl=TRUE)
  out_txt <- paste0(base, ".annotated.txt")
} else {
  out_txt <- opt$options$output
}

# ---- Files ----
intron_bed     <- paste0(anno_prefix, "_all_introns.bed.gz")
fiveprime_bed  <- paste0(anno_prefix, "_fiveprime.bed.gz")
threeprime_bed <- paste0(anno_prefix, "_threeprime.bed.gz")

for (f in c(counts_gz, intron_bed, fiveprime_bed, threeprime_bed)) {
  if (!file.exists(f)) stop("Missing file: ", f)
}

# ---- Helpers ----
add_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

# Verdict logic (per intron), mirroring LeafViz
verdict_for_intron <- function(row_one, gene_strand, fprime_tbl, tprime_tbl, bothSS_tbl) {
  chr <- row_one$chr; st <- row_one$start; en <- row_one$end

  fprime_intron <- fprime_tbl %>% dplyr::filter(chr == !!chr, start == !!st)
  tprime_intron <- tprime_tbl %>% dplyr::filter(chr == !!chr, end   == !!en)
  bothSS_intron <- bothSS_tbl %>% dplyr::filter(chr == !!chr, start == !!st, end == !!en)

  unknown_3p <- nrow(tprime_intron) == 0 || all(is.na(tprime_intron$gene))
  unknown_5p <- nrow(fprime_intron) == 0 || all(is.na(fprime_intron$gene))

  if (is.na(gene_strand)) return("unknown_strand")
  if (unknown_3p && unknown_5p) return("cryptic_unanchored")

  if ((unknown_3p && !unknown_5p && gene_strand == "+") ||
      (unknown_5p && !unknown_3p && gene_strand == "-")) return("cryptic_threeprime")

  if ((!unknown_3p && unknown_5p && gene_strand == "+") ||
      (!unknown_5p && unknown_3p && gene_strand == "-")) return("cryptic_fiveprime")

  both_sites_annot <- (!unknown_3p) && (!unknown_5p)
  if (both_sites_annot) {
    if (nrow(bothSS_intron) > 0 && all(!is.na(bothSS_intron$gene))) {
      return("annotated")
    } else {
      return("novel annotated pair")
    }
  }
  "error"
}

# ---- Load annotation DBs ----
message("Loading annotation data...")
intron_db <- fread(cmd = paste("zcat < ", shQuote(intron_bed)),
                   sep = "\t", header = FALSE, data.table = FALSE)
colnames(intron_db)[1:4] <- c("chr","start","end","gene")
intron_db$chr <- add_chr(intron_db$chr)

fiveprime_db <- fread(cmd = paste("zcat < ", shQuote(fiveprime_bed)),
                      sep = "\t", header = FALSE, data.table = FALSE)
colnames(fiveprime_db)[1:7] <- c("chr","start","end","gene","gene_id","strand","transcript")
fiveprime_db$chr <- add_chr(fiveprime_db$chr)

threeprime_db <- fread(cmd = paste("zcat < ", shQuote(threeprime_bed)),
                       sep = "\t", header = FALSE, data.table = FALSE)
colnames(threeprime_db)[1:7] <- c("chr","start","end","gene","gene_id","strand","transcript")
threeprime_db$chr <- add_chr(threeprime_db$chr)

# ---- Parse all junctions from counts file (TAB-separated) ----
message("Parsing junction IDs from counts...")

# Read only column 1 (junction ID); header row comes first and is removed
tmp <- fread(cmd = paste("zcat < ", shQuote(counts_gz)),
             header = FALSE, sep = "\t", select = 1,
             colClasses = "character", data.table = FALSE)

if (nrow(tmp) < 2) stop("Counts file has no data rows.")

# First row is sample header (starts with a leading tab) -> drop it
tmp <- tmp[-1, , drop = FALSE]

# First column is the junction ID
junction_ids <- trimws(tmp[[1]])
# Safety: ensure just the first token
junction_ids <- sub("[ \t].*$", "", junction_ids, perl = TRUE)

# Build junction table
split4 <- str_split_fixed(junction_ids, ":", 4)
junctions <- data.frame(
  intron_id = junction_ids,
  chr       = add_chr(split4[,1]),
  start     = as.numeric(split4[,2]),
  end       = as.numeric(split4[,3]),
  clusterID = split4[,4],
  stringsAsFactors = FALSE
)

# ---- Precompute per-junction annotation ("verdict") ----
message("Annotating junctions... (this mirrors your LeafViz logic)")

# Lookups
bothSS_all       <- intron_db %>% dplyr::select(chr, start, end, gene)
fiveprime_keyed  <- fiveprime_db %>% dplyr::select(chr, start, gene, gene_id, strand, transcript)
threeprime_keyed <- threeprime_db %>% dplyr::select(chr, start, gene, gene_id, strand, transcript)

verdict_vec <- character(nrow(junctions))
names(verdict_vec) <- junctions$intron_id

clusters <- unique(junctions$clusterID)
pb_total <- length(clusters); i <- 0L

for (clu in clusters) {
  i <- i + 1L
  if (i %% 1000 == 0) message(sprintf("  ... %d / %d clusters", i, pb_total))

  cluster_df <- junctions %>% dplyr::filter(clusterID == !!clu)

  # Per-cluster tables
  fprimeClu <- cluster_df %>%
    dplyr::left_join(fiveprime_keyed, by = c("chr","start"))
  tprimeClu <- cluster_df %>%
    dplyr::left_join(threeprime_keyed, by = c("chr"="chr", "end"="start")) # end ↔ 3' start
  bothSSClu <- cluster_df %>%
    dplyr::left_join(bothSS_all, by = c("chr","start","end"))

  # Cluster gene and strand consensus
  cluster_gene <- names(sort(table(c(tprimeClu$gene, fprimeClu$gene)), decreasing = TRUE))[1]
  if (length(cluster_gene) == 0 || is.na(cluster_gene)) cluster_gene <- "."

  gene_strand <- NA
  if (cluster_gene != ".") {
    strands <- c(tprimeClu$strand, fprimeClu$strand)
    strands <- unique(strands[!is.na(strands) & strands != "."])
    if (length(strands) == 1) gene_strand <- strands[1]
  }

  # Verdict per intron in this cluster
  for (r in seq_len(nrow(cluster_df))) {
    intron_row <- cluster_df[r,]
    v <- verdict_for_intron(intron_row, gene_strand, fprimeClu, tprimeClu, bothSSClu)
    verdict_vec[intron_row$intron_id] <- v
  }
}

# Any missing verdicts -> "."
verdict_vec[is.na(verdict_vec) | verdict_vec == ""] <- "."

# Optional sanity check: parsed junction IDs must equal names(verdict_vec)
if (!setequal(names(verdict_vec), junction_ids)) {
  missing_in_verdict <- setdiff(junction_ids, names(verdict_vec))
  extra_in_verdict   <- setdiff(names(verdict_vec), junction_ids)
  warning(sprintf("Mismatch between junction sets. Missing: %d, Extra: %d",
                  length(missing_in_verdict), length(extra_in_verdict)))
}

# ---- Stream input and write output with trailing 'annotation' (TAB) ----
message("Writing annotated counts to: ", out_txt)

con_in  <- gzfile(counts_gz, open = "rt"); on.exit(try(close(con_in),  silent=TRUE), add = TRUE)
con_out <- file(out_txt,       open = "wt"); on.exit(try(close(con_out), silent=TRUE), add = TRUE)

# Header: original + TAB + annotation
hdr <- readLines(con_in, n = 1L)
if (length(hdr) == 0) stop("Empty counts file.")
writeLines(paste0(hdr, "\tannotation"), con_out)

repeat {
  lines <- readLines(con_in, n = 10000L)
  if (!length(lines)) break

  # First token = junction id (trim leading whitespace)
  ids <- sub("^([^\t ]+).*", "\\1", trimws(lines))
  ann <- verdict_vec[ids]                 # name-based lookup; missing -> NA
  ann[is.na(ann) | !nzchar(ann)] <- "."  # replace unknown/blank with '.'

  # Append with TAB, preserve TSV
  writeLines(paste0(lines, "\t", ann), con_out)
}

message("Done.")
