#!/usr/bin/env Rscript
# 02_de.R — DESeq2 on the six-group design, two contrasts, shrunken LFCs.
#
# Design is ~group (cell x stage interaction folded into one factor), which
# lets any pair of the six groups be contrasted directly. Two are reported:
#   B_vs_L_virgin : basal vs luminal within virgin animals (cell-type axis)
#   lact_vs_preg_luminal : lactating vs pregnant luminal cells (stage axis)
#
# LFCs are shrunk with apeglm; the unshrunken Wald statistics are kept for
# ranking in the enrichment step.

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(readr)
})

RES <- "results"
ALPHA <- 0.05
LFC_MIN <- 1          # |log2FC| > 1 for the "strong" call, reported alongside alpha

run_de <- function(dds) {
  dds <- DESeq(dds, quiet = TRUE)
  dds
}

contrast_table <- function(dds, num, den, name) {
  res <- results(dds, contrast = c("group", num, den), alpha = ALPHA)
  # apeglm needs a coefficient, so relevel and refit the shrinkage on the pair
  dds2 <- dds
  colData(dds2)$group <- relevel(droplevels(factor(colData(dds2)$group)), ref = den)
  dds2 <- nbinomWaldTest(estimateDispersions(
    DESeqDataSetFromMatrix(counts(dds2), colData(dds2), ~ group), quiet = TRUE), quiet = TRUE)
  coefname <- paste0("group_", num, "_vs_", den)
  shr <- lfcShrink(dds2, coef = coefname, type = "apeglm", quiet = TRUE)

  out <- as.data.frame(res) %>%
    tibble::rownames_to_column("entrez") %>%
    mutate(lfc_shrunk = as.data.frame(shr)[entrez, "log2FoldChange"],
           contrast = name) %>%
    arrange(padj)
  out
}

summarise_contrast <- function(tab) {
  # Significance comes from the unshrunken Wald test; the effect-size call uses
  # the SHRUNKEN LFC, because that is the axis reported in the volcano plots.
  # The two disagree by 150-300 genes per contrast (shrinkage moves
  # low-information genes across the |LFC|=1 line), so the pipeline must pick
  # one and use it everywhere. Both counts are returned to make that explicit.
  sig <- !is.na(tab$padj) & tab$padj < ALPHA
  strong_shrunk   <- sig & abs(tab$lfc_shrunk)     > LFC_MIN
  strong_unshrunk <- sig & abs(tab$log2FoldChange) > LFC_MIN
  list(n_tested = sum(!is.na(tab$padj)),
       n_sig = sum(sig),
       n_strong = sum(strong_shrunk),
       n_strong_unshrunk = sum(strong_unshrunk),
       n_up = sum(strong_shrunk & tab$lfc_shrunk > 0),
       n_down = sum(strong_shrunk & tab$lfc_shrunk < 0))
}

# Biological positive control for the lactating-vs-pregnant luminal contrast.
#
# The sample-sheet tests already catch a positional join. They do NOT catch a
# swapped contrast: contrast_table(dds, "luminal_pregnant", "luminal_lactating")
# is a valid DESeq2 call that inverts every sign, and every count in
# de_summary.csv is sign-agnostic (n_sig, n_strong) or symmetric under the swap
# (n_up and n_down exchange), so nothing downstream notices.
#
# What pins the direction is biology outside the count matrix: milk-protein
# genes are induced in lactation. Asserting that they come out UP turns an
# inverted comparison into a failed run rather than a published figure with
# every arrow backwards.
#
# The markers are named by symbol and resolved through the annotation, so the
# check cannot be satisfied by whichever genes happen to top the table.
POS_CONTROL_MARKERS <- c("Csn2", "Wap", "Lalba", "Csn1s1", "Csn1s2a")

positive_control <- function(tab, sym, markers = POS_CONTROL_MARKERS,
                             min_lfc = 1, max_padj = 1e-4) {
  ids <- names(sym)[match(markers, unname(sym))]
  missing <- markers[is.na(ids)]
  if (length(missing))
    stop("positive control: marker not in annotation: ", paste(missing, collapse = ", "))

  i <- match(ids, tab$entrez)
  if (anyNA(i))
    stop("positive control: marker not in the tested set: ",
         paste(markers[is.na(i)], collapse = ", "))

  out <- data.frame(symbol = markers, entrez = ids,
                    lfc_shrunk = tab$lfc_shrunk[i], padj = tab$padj[i],
                    stringsAsFactors = FALSE)
  bad <- is.na(out$lfc_shrunk) | out$lfc_shrunk < min_lfc |
         is.na(out$padj) | out$padj > max_padj
  if (any(bad))
    stop("positive control FAILED — lactation markers are not up in this ",
         "contrast, so the comparison is probably inverted:\n",
         paste(sprintf("  %-8s lfc_shrunk=%+.2f padj=%.3g",
                       out$symbol[bad], out$lfc_shrunk[bad], out$padj[bad]),
               collapse = "\n"))
  out
}

if (sys.nframe() == 0) {
  dds <- readRDS(file.path(RES, "dds_raw.rds"))
  dds <- run_de(dds)
  saveRDS(dds, file.path(RES, "dds_fitted.rds"))

  contrasts <- list(
    B_vs_L_virgin        = c("basal_virgin", "luminal_virgin"),
    lact_vs_preg_luminal = c("luminal_lactating", "luminal_pregnant")
  )
  all_tabs <- lapply(names(contrasts), function(nm) {
    p <- contrasts[[nm]]
    contrast_table(dds, p[1], p[2], nm)
  })
  names(all_tabs) <- names(contrasts)
  saveRDS(all_tabs, file.path(RES, "de_tables.rds"))

  for (nm in names(all_tabs)) {
    s <- summarise_contrast(all_tabs[[nm]])
    cat(sprintf("%-22s tested=%d sig(padj<%.2f)=%d strong(|LFC|>%d)=%d up=%d down=%d\n",
                nm, s$n_tested, ALPHA, s$n_sig, LFC_MIN, s$n_strong, s$n_up, s$n_down))
    write_csv(all_tabs[[nm]], file.path(RES, paste0("de_", nm, ".csv")))
  }
}
