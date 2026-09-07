library(testthat)
suppressPackageStartupMessages({library(DESeq2); library(jsonlite)})

# testthat runs test files with the working directory set to tests/, so
# resolve the project root explicitly and read everything relative to it.
ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".."))
if (!dir.exists(file.path(ROOT, "R")))
  stop("cannot locate project root; run: Rscript tests/run_tests.R")
setwd(ROOT)

source("R/01_load.R")
source("R/00_annotation.R")
source("R/02_de.R")
source("R/06_sensitivity.R")

oracle <- fromJSON("tests/oracle.json")

test_that("counts matrix has the published shape", {
  cd <- load_counts()
  expect_equal(nrow(cd$counts), oracle$n_genes_raw)
  expect_equal(ncol(cd$counts), oracle$n_samples)
  expect_true(all(grepl("^MCL1-[DL][A-L]$", colnames(cd$counts))))
})

# THE regression test for this dataset. The counts file lists basal libraries
# first (MCL1-DG..DL); the GEO series matrix lists luminal first. Joining by
# position instead of by label silently swaps the two cell populations and
# inverts every fold change, with no error anywhere downstream.
test_that("sample sheet joins by label, not by position", {
  cd <- load_counts()
  ss <- load_sample_sheet(cd$counts)
  expect_identical(as.character(ss$label), colnames(cd$counts))
  expect_true(all(grepl("^MCL1-D", ss$label[ss$cell == "basal"])))
  expect_true(all(grepl("^MCL1-L", ss$label[ss$cell == "luminal"])))
  # a positionally-joined sheet would fail this: the first counts column is
  # basal, but the first series-matrix sample is luminal
  expect_equal(as.character(ss$cell[1]), "basal")
})

test_that("a positionally-joined sheet is rejected", {
  cd <- load_counts()
  ss <- read.csv("data/sample_sheet.csv")
  bad <- ss
  bad$label <- rev(bad$label)          # scramble the link
  f <- tempfile(fileext = ".csv"); write.csv(bad, f, row.names = FALSE)
  expect_error(load_sample_sheet(cd$counts, f), regexp = NULL)
  unlink(f)
})

test_that("gene filter is reproducible", {
  cd <- load_counts()
  expect_equal(nrow(filter_genes(cd$counts)), oracle$n_genes_kept)
})

test_that("library sizes match the independent Python oracle", {
  cd <- load_counts()
  ls <- colSums(cd$counts)
  expect_equal(unname(max(ls)), oracle$libsize_max)
  expect_equal(unname(min(ls)), oracle$libsize_min)
  expect_equal(names(which.max(ls)), oracle$libsize_max_sample)
})

test_that("design is balanced 2 per group", {
  cd <- load_counts()
  ss <- load_sample_sheet(cd$counts)
  expect_true(all(table(ss$group) == 2))
  expect_equal(nlevels(droplevels(ss$group)), 6)
})

test_that("annotation covers almost all tested genes", {
  skip_if_not(file.exists("data/annotation/Mus_musculus.gene_info.gz"))
  sym <- load_gene_symbols()
  cd <- load_counts()
  kept <- rownames(filter_genes(cd$counts))
  expect_gt(mean(kept %in% names(sym)), 0.95)
})

test_that("milk-gene exclusion set is defined by family, not by expression", {
  skip_if_not(file.exists("data/annotation/Mus_musculus.gene_info.gz"))
  sym <- load_gene_symbols()
  cd <- load_counts()
  milk <- milk_gene_ids(sym, rownames(cd$counts))
  expect_true(length(milk) >= 15)
  expect_true(all(c("Csn2", "Wap", "Lalba") %in% unname(sym[milk])))
  # the set must not depend on which libraries are in the object
  expect_identical(milk, milk_gene_ids(sym, rownames(cd$counts)))
})

# The contrast-direction tests below fit DESeq2 (~45 s on this design). They are
# the only tests that touch the model, and they are the ones that catch a
# swapped contrast, which no count-based check can see.
test_that("lactation markers are up, and a swapped contrast fails", {
  skip_if_not(file.exists("data/annotation/Mus_musculus.gene_info.gz"))
  sym <- load_gene_symbols()
  dds <- run_de(build_dds())

  fwd <- contrast_table(dds, "luminal_lactating", "luminal_pregnant", "fwd")
  pc <- positive_control(fwd, sym)
  expect_equal(nrow(pc), length(POS_CONTROL_MARKERS))
  expect_true(all(pc$lfc_shrunk > 1))
  expect_true(all(pc$padj < 1e-4))

  # falsification: the same call with numerator and denominator exchanged is
  # valid DESeq2 and inverts every sign, so the gate must reject it
  rev <- contrast_table(dds, "luminal_pregnant", "luminal_lactating", "rev")
  expect_error(positive_control(rev, sym), "positive control FAILED")

  # and the markers must be named, not discovered: a marker absent from the
  # annotation is an error, never a silent pass
  expect_error(positive_control(fwd, sym, markers = c("Csn2", "NotAGene")),
               "not in annotation")
})
