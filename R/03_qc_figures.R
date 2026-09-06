#!/usr/bin/env Rscript
# 03_qc_figures.R — QC panels that justify the normalisation choice.
#
# The point of the first panel: in the two lactating luminal libraries, ~54%
# of all reads sit in just ten genes — eight of them caseins/whey proteins,
# plus COX1 and Scd1. Their raw library sizes are the LARGEST of the twelve,
# yet their DESeq2 size factors are the SMALLEST. Normalising by total counts
# would therefore deflate every other gene in those samples. This is the
# concrete reason median-of-ratios is used instead, and it is measurable in
# this dataset rather than asserted.

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(dplyr); library(tidyr)
})

FIG <- "results/figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

BLUE <- "#1F5FA8"; RED <- "#8C2D2D"; GREY <- "#5A5A5A"

qc_table <- function(dds) {
  cm <- counts(dds)
  data.frame(
    label     = colnames(cm),
    cell      = colData(dds)$cell,
    stage     = colData(dds)$stage,
    libsize   = colSums(cm),
    sizefac   = sizeFactors(dds),
    top10     = apply(cm, 2, function(x) sum(sort(x, decreasing = TRUE)[1:10]) / sum(x))
  )
}

# top_syms: symbols of the top-10 genes in the most skewed library, and which
# of them are milk-protein family. Passed in rather than hardcoded so the
# subtitle cannot drift from the data — an earlier version claimed all ten
# were milk-protein genes when two (COX1, Scd1) are not.
plot_libsize_vs_sizefactor <- function(qc, top_syms = NULL, top_is_milk = NULL,
                                       path = file.path(FIG, "qc_normalisation.png")) {
  # label the samples that carry the argument: the two most skewed and,
  # for contrast, the least skewed one
  qc$flag <- rank(-qc$top10) <= 2 | rank(qc$top10) == 1
  lab <- qc[qc$flag, ]
  lab$txt <- sprintf("%s\n%s %s\n%.0f%% in top 10", lab$label, lab$cell,
                     ifelse(lab$stage == "2 day lactation", "lactating",
                            as.character(lab$stage)), 100 * lab$top10)

  # every number and gene name in the subtitle is measured here, not asserted
  worst <- qc[order(-qc$top10), ][1:2, ]
  sub <- if (is.null(top_syms)) {
    sprintf("Two lactating luminal samples hold %.0f%% and %.0f%% of their reads in ten genes;\nscaling by total counts would deflate every other gene in them, so median-of-ratios\nis used instead",
            100 * worst$top10[1], 100 * worst$top10[2])
  } else {
    sprintf("%s holds %.0f%% of its reads in ten genes — %d caseins/whey (%s) plus %s;\nscaling by total counts would deflate every other gene in it, so median-of-ratios\nis used instead",
            worst$label[1], 100 * worst$top10[1], sum(top_is_milk),
            paste(head(top_syms[top_is_milk], 3), collapse = ", "),
            paste(top_syms[!top_is_milk], collapse = " and "))
  }

  p <- ggplot(qc, aes(libsize / 1e6, sizefac)) +
    geom_point(aes(colour = top10, shape = cell), size = 3) +
    ggrepel::geom_text_repel(data = lab, aes(label = txt), size = 2.4, colour = GREY,
                             lineheight = 0.95, min.segment.length = 0,
                             box.padding = 0.7, point.padding = 0.4,
                             segment.colour = GREY, segment.size = 0.3, seed = 1) +
    scale_colour_gradient(low = BLUE, high = RED, name = "share of reads\nin top 10 genes",
                          labels = scales::percent_format(accuracy = 1)) +
    scale_shape_manual(values = c(luminal = 16, basal = 17), name = "cell type") +
    scale_y_continuous(expand = expansion(mult = 0.14)) +
    scale_x_continuous(expand = expansion(mult = 0.10)) +
    labs(x = "Raw library size (million reads)",
         y = "DESeq2 size factor",
         title = "The largest libraries get the smallest size factors",
         subtitle = sub) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(hjust = 0, face = "plain"),
          plot.subtitle = element_text(hjust = 0, colour = GREY, size = 7))
  ggsave(path, plot = p, width = 6.2, height = 4.2, dpi = 300)
  path
}

plot_pca <- function(dds, path = file.path(FIG, "qc_pca.png")) {
  vsd <- vst(dds, blind = TRUE)
  pc <- plotPCA(vsd, intgroup = c("cell", "stage"), returnData = TRUE)
  pv <- round(100 * attr(pc, "percentVar"))
  p <- ggplot(pc, aes(PC1, PC2, colour = stage, shape = cell)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(luminal = 16, basal = 17)) +
    labs(x = sprintf("PC1 (%d%% variance)", pv[1]),
         y = sprintf("PC2 (%d%% variance)", pv[2]),
         title = "Cell type separates on PC1, developmental stage on PC2",
         colour = "stage", shape = "cell type") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(hjust = 0, face = "plain"))
  ggsave(path, plot = p, width = 5.4, height = 3.6, dpi = 300)
  list(path = path, percentVar = pv, data = pc)
}

if (sys.nframe() == 0) {
  dds <- readRDS("results/dds_fitted.rds")
  qc <- qc_table(dds)
  write.csv(qc, "results/qc_table.csv", row.names = FALSE)
  plot_libsize_vs_sizefactor(qc)
  pca <- plot_pca(dds)
  print(qc[order(-qc$top10), ], digits = 3)
  cat("PC variance:", pca$percentVar, "\n")
}
