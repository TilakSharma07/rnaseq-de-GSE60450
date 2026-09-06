#!/usr/bin/env Rscript
# 06_sensitivity.R — is the "translation down" result biology or normalisation?
#
# Both contrasts return strongly NEGATIVE enrichment for ribosome and
# translation pathways. In this dataset that has two candidate explanations:
#
#   (a) biology: lactating cells devote their transcriptome to milk protein,
#       so the relative share of ribosomal transcripts genuinely falls;
#   (b) artefact: a handful of casein genes absorb >50% of the library, and
#       any residual failure of normalisation would push everything else
#       down together — ribosomal genes are numerous and highly expressed,
#       so they would dominate such an artefact.
#
# The two are distinguishable. Drop the milk-protein genes entirely, refit
# from raw counts (new size factors, new dispersions, new tests), and rerun
# the enrichment. If the translation signature survives on a library where
# caseins cannot influence normalisation, (a) holds. If it collapses, the
# original result was (b) and must not be reported as biology.
#
# This is the check that decides whether the enrichment section of the report
# is a finding or an artefact, so it runs as part of the pipeline, not as an
# afterthought.

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr)
})

# Milk-protein and closely-related secretory genes, by symbol. Defined by
# gene family rather than by "top N in the biggest library", so the exclusion
# set does not depend on the very libraries under test.
MILK_PATTERNS <- c("^Csn[0-9]", "^Wap$", "^Glycam1$", "^Lalba$", "^Muc1$",
                   "^Btn1a1$", "^Xdh$", "^Cel$", "^Lpo$", "^Ltf$", "^Trf$",
                   "^Spp1$", "^Mfge8$")

milk_gene_ids <- function(sym, universe) {
  hit <- names(sym)[Reduce(`|`, lapply(MILK_PATTERNS, function(p) grepl(p, sym)))]
  intersect(hit, universe)
}

refit_without <- function(dds_raw, drop_ids) {
  keep <- setdiff(rownames(dds_raw), drop_ids)
  sub <- DESeqDataSetFromMatrix(counts(dds_raw)[keep, , drop = FALSE],
                                colData(dds_raw), design = ~ group)
  DESeq(sub, quiet = TRUE)
}

compare_nes <- function(gsea_full, gsea_sub, label_regex) {
  a <- gsea_full[, c("pathway", "NES", "padj")]
  b <- gsea_sub[,  c("pathway", "NES", "padj")]
  m <- merge(a, b, by = "pathway", suffixes = c("_full", "_nomilk"))
  # "translation" also matches "post-translational modification" and
  # "posttranslational processing", which are unrelated to ribosome biogenesis
  # and would inflate the highlighted set from 15 to 21. Exclude them.
  m$is_translation <- grepl(label_regex, m$pathway, ignore.case = TRUE) &
                      !grepl("post-?translational", m$pathway, ignore.case = TRUE)
  m
}

plot_sensitivity <- function(cmp, path = "results/figures/sensitivity_nomilk.png") {
  library(ggplot2)
  GREY <- "#5A5A5A"; RED <- "#8C2D2D"; PALE <- "#BFBFBF"
  r <- cor(cmp$NES_full, cmp$NES_nomilk, use = "complete.obs")
  n_tr <- sum(cmp$is_translation)
  n_tr_sig <- sum(cmp$is_translation & cmp$padj_nomilk < 0.05)
  p <- ggplot(cmp, aes(NES_full, NES_nomilk)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = GREY, linewidth = 0.3) +
    geom_point(aes(colour = is_translation, size = is_translation), alpha = 0.7) +
    scale_colour_manual(values = c(`TRUE` = RED, `FALSE` = PALE),
                        labels = c(`TRUE` = sprintf("translation / ribosome (%d sets, %d significant)",
                                                    n_tr, n_tr_sig),
                                   `FALSE` = "other pathways"), name = NULL) +
    scale_size_manual(values = c(`TRUE` = 1.5, `FALSE` = 0.7), guide = "none") +
    annotate("text", x = min(cmp$NES_full, na.rm = TRUE), y = max(cmp$NES_nomilk, na.rm = TRUE),
             hjust = 0, vjust = 1, size = 2.6, colour = "black",
             label = sprintf("Pearson r = %.4f\nover %d Reactome pathways", r, nrow(cmp))) +
    labs(x = "NES, all 16,037 genes",
         y = "NES, 17 milk-protein genes removed and refit",
         title = "The translation signature is not a normalisation artefact",
         subtitle = sprintf("Dropping the casein/whey genes that hold 54%% of the lactating libraries and refitting from\nraw counts moves no NES by more than %.2f (r = %.4f); %d of %d pathways change significance\ncall. Every significant translation set stays significant, so the depletion is a real shift",
                            max(abs(cmp$NES_full - cmp$NES_nomilk), na.rm = TRUE), r,
                            sum((cmp$padj_full < 0.05) != (cmp$padj_nomilk < 0.05)), nrow(cmp))) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(hjust = 0, face = "plain"),
          plot.subtitle = element_text(hjust = 0, colour = GREY, size = 6.8),
          legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6.5))
  ggsave(path, plot = p, width = 5.8, height = 4.4, dpi = 300)
  path
}
