#!/usr/bin/env Rscript
# 00_annotation.R — Entrez -> symbol, and Entrez -> Reactome pathway sets.
#
# Deliberately does NOT use org.Mm.eg.db. That package is a versioned snapshot
# whose contents change between Bioconductor releases, so a pipeline pinned to
# it is only reproducible if the exact release is also pinned. Here the two
# source files are downloaded from NCBI and Reactome, their release identity
# is recorded, and both are small enough to cache next to the analysis.
#
#   Mus_musculus.gene_info.gz  NCBI Gene, authoritative Entrez -> symbol
#   NCBI2Reactome.txt          Reactome, gene -> pathway with evidence codes
#
# Reactome rows carry an evidence code; IEA (electronic annotation) rows are
# kept separately from TAS (traceable author statement) so an enrichment run
# can choose. Default keeps both and records the split.

suppressPackageStartupMessages({
  library(readr); library(dplyr)
})

ANN <- "data/annotation"

load_gene_symbols <- function(path = file.path(ANN, "Mus_musculus.gene_info.gz")) {
  gi <- read_tsv(path, show_col_types = FALSE, progress = FALSE,
                 col_types = cols(.default = col_character()))
  names(gi)[1] <- "tax_id"
  stopifnot(all(c("GeneID", "Symbol", "type_of_gene") %in% names(gi)))
  gi <- gi[gi$tax_id == "10090", ]           # mouse only, in case the file changes
  setNames(gi$Symbol, gi$GeneID)
}

load_gene_descriptions <- function(path = file.path(ANN, "Mus_musculus.gene_info.gz")) {
  gi <- read_tsv(path, show_col_types = FALSE, progress = FALSE,
                 col_types = cols(.default = col_character()))
  names(gi)[1] <- "tax_id"
  gi <- gi[gi$tax_id == "10090", ]
  setNames(gi$description, gi$GeneID)
}

# Reactome pathway sets, restricted to a gene universe (the tested genes)
load_reactome <- function(universe, path = file.path(ANN, "NCBI2Reactome.txt"),
                          species = "Mus musculus", min_size = 10, max_size = 500) {
  rc <- read_tsv(path, col_names = c("gene", "pathway_id", "url", "pathway", "evidence", "species"),
                 show_col_types = FALSE, progress = FALSE,
                 col_types = cols(.default = col_character()))
  rc <- rc[rc$species == species, ]
  rc$pathway <- trimws(rc$pathway)
  ev <- table(rc$evidence)
  message(sprintf("reactome: %d %s rows, evidence: %s",
                  nrow(rc), species,
                  paste(sprintf("%s=%d", names(ev), as.integer(ev)), collapse = " ")))

  rc <- rc[rc$gene %in% universe, ]
  sets <- split(rc$gene, rc$pathway)
  sets <- lapply(sets, unique)
  sizes <- lengths(sets)
  sets <- sets[sizes >= min_size & sizes <= max_size]
  message(sprintf("reactome: %d pathways with %d-%d genes inside the tested universe",
                  length(sets), min_size, max_size))
  sets
}

annotation_provenance <- function() {
  f1 <- file.path(ANN, "Mus_musculus.gene_info.gz")
  f2 <- file.path(ANN, "NCBI2Reactome.txt")
  data.frame(
    file = c(basename(f1), basename(f2)),
    source = c("NCBI Gene FTP", "Reactome current release"),
    bytes = c(file.size(f1), file.size(f2)),
    downloaded = format(c(file.mtime(f1), file.mtime(f2)), "%Y-%m-%d")
  )
}
