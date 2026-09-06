#!/usr/bin/env Rscript
# build_sample_sheet.R — derive data/sample_sheet.csv from the GEO series matrix.
#
# Parses !Sample_title, !Sample_geo_accession, !Sample_description (which holds
# "Sample name: MCL1-xx", the only field linking GEO samples to counts columns)
# and the two !Sample_characteristics_ch1 rows carrying immunophenotype and
# developmental stage. Emits one row per counts column, IN counts-column order.

suppressPackageStartupMessages(library(readr))

SM <- "data/raw/GSE60450_series_matrix.txt.gz"
CT <- "data/raw/GSE60450_counts.txt.gz"
OUT <- "data/sample_sheet.csv"

lines <- readLines(gzfile(SM), warn = FALSE)
field <- function(prefix) {
  hit <- lines[startsWith(lines, prefix)]
  if (!length(hit)) stop("field not found in series matrix: ", prefix)
  lapply(hit, function(L) gsub('"', '', strsplit(L, "\t")[[1]][-1]))
}

gsm   <- field("!Sample_geo_accession")[[1]]
title <- field("!Sample_title")[[1]]
label <- sub("^Sample name: ", "", field("!Sample_description")[[1]])

# both immunophenotype and developmental stage arrive as !Sample_characteristics_ch1
chars <- field("!Sample_characteristics_ch1")
pick <- function(key) {
  for (v in chars) if (startsWith(v[1], key)) return(sub(paste0("^", key, ": "), "", v))
  stop("characteristic not found: ", key)
}
cell  <- sub(" cell population$", "", pick("immunophenotype"))
stage <- pick("developmental stage")

hdr  <- strsplit(readLines(gzfile(CT), n = 1), "\t")[[1]]
cols <- sub("_.*$", "", hdr[-(1:2)])

stopifnot(length(gsm) == length(cols), setequal(label, cols))
ord <- match(cols, label)          # join by label; NEVER by position
ss <- data.frame(
  gsm = gsm[ord], label = label[ord], title = title[ord],
  cell = cell[ord], stage = stage[ord],
  group = paste0(cell[ord], "_", c(virgin = "virgin",
                                   `18.5 day pregnancy` = "pregnant",
                                   `2 day lactation` = "lactating")[stage[ord]]))
stopifnot(identical(ss$label, cols),
          all(grepl("^MCL1-D", ss$label[ss$cell == "basal"])),
          all(grepl("^MCL1-L", ss$label[ss$cell == "luminal"])),
          all(table(ss$group) == 2))
write_csv(ss, OUT)
cat(sprintf("  wrote %s: %d samples, %d groups\n", OUT, nrow(ss), length(unique(ss$group))))
print(table(ss$cell, ss$stage))
