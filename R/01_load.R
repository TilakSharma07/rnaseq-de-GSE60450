#!/usr/bin/env Rscript
# 01_load.R — read GSE60450 counts, attach the sample sheet, filter low counts.
#
# The one trap in this dataset: the counts matrix and the GEO series matrix
# list the 12 samples in DIFFERENT orders. The counts file starts with the
# basal libraries (MCL1-DG..DL); the series matrix starts with the luminal
# ones (GSM1480291 = MCL1-LA, luminal). Assigning sample titles by position
# would silently swap basal and luminal for all 12 samples and invert every
# log2 fold change. So we join on the library label and assert the result.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
})

DATA_RAW <- "data/raw"
DATA     <- "data"
RES      <- "results"
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

load_counts <- function(path = file.path(DATA_RAW, "GSE60450_counts.txt.gz")) {
  raw <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  stopifnot("EntrezGeneID" %in% names(raw), "Length" %in% names(raw))
  gene_len <- setNames(raw$Length, as.character(raw$EntrezGeneID))
  mat <- as.matrix(raw[, setdiff(names(raw), c("EntrezGeneID", "Length"))])
  rownames(mat) <- as.character(raw$EntrezGeneID)
  # column headers are "MCL1-DG_BC2CTUACXX_ACTTGA_L002_R1" -> keep the label
  colnames(mat) <- sub("_.*$", "", colnames(mat))
  list(counts = mat, gene_length = gene_len)
}

load_sample_sheet <- function(counts, path = file.path(DATA, "sample_sheet.csv")) {
  ss <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  stopifnot(all(c("gsm", "label", "cell", "stage", "group") %in% names(ss)))

  # join on label, never on position
  if (!setequal(ss$label, colnames(counts)))
    stop("sample sheet labels do not match counts columns")
  ss <- ss[match(colnames(counts), ss$label), ]

  # the assertion that makes the ordering bug impossible to reintroduce:
  # every basal library must be one of DG..DL, every luminal one LA..LF.
  basal_ids   <- ss$label[ss$cell == "basal"]
  luminal_ids <- ss$label[ss$cell == "luminal"]
  stopifnot(all(grepl("^MCL1-D", basal_ids)), all(grepl("^MCL1-L", luminal_ids)))
  stopifnot(identical(ss$label, colnames(counts)))

  ss$cell  <- factor(ss$cell,  levels = c("luminal", "basal"))
  ss$stage <- factor(ss$stage, levels = c("virgin", "18.5 day pregnancy", "2 day lactation"))
  ss$group <- factor(ss$group)
  ss
}

# keep genes with >= MIN_COUNT reads in at least MIN_SAMPLES libraries
MIN_COUNT   <- 10
MIN_SAMPLES <- 2

filter_genes <- function(counts, min_count = MIN_COUNT, min_samples = MIN_SAMPLES) {
  keep <- rowSums(counts >= min_count) >= min_samples
  message(sprintf("filter: %d of %d genes kept (>= %d reads in >= %d samples)",
                  sum(keep), nrow(counts), min_count, min_samples))
  counts[keep, , drop = FALSE]
}

build_dds <- function() {
  cd <- load_counts()
  ss <- load_sample_sheet(cd$counts)
  mat <- filter_genes(cd$counts)
  dds <- DESeqDataSetFromMatrix(countData = mat, colData = ss, design = ~ group)
  metadata(dds)$gene_length <- cd$gene_length[rownames(mat)]
  dds
}

if (sys.nframe() == 0) {
  dds <- build_dds()
  saveRDS(dds, file.path(RES, "dds_raw.rds"))
  cat(sprintf("dds: %d genes x %d samples\n", nrow(dds), ncol(dds)))
  print(table(colData(dds)$cell, colData(dds)$stage))
}
