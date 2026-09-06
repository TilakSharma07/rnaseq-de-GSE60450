#!/usr/bin/env Rscript
# 04_de_figures.R — volcano and MA panels per contrast, gene-symbol labelled.
#
# Volcano y-axis uses the unshrunken p-value (significance) against the
# SHRUNKEN log2 fold change (effect size). Mixing them is deliberate: apeglm
# shrinkage moves low-information genes toward zero without touching their
# evidence, so the shrunken axis is the honest effect axis while the p-value
# stays the test result. Labelled genes are the top by padj, not by fold
# change, so a large-but-noisy gene never gets named.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ggrepel)
})

FIG <- "results/figures"
ALPHA <- 0.05; LFC_MIN <- 1
BLUE <- "#1F5FA8"; RED <- "#8C2D2D"; GREY <- "#5A5A5A"; PALE <- "#BFBFBF"

label_genes <- function(tab, sym, n = 12) {
  t <- tab[!is.na(tab$padj) & tab$padj < ALPHA & abs(tab$lfc_shrunk) > LFC_MIN, ]
  t <- t[order(t$padj), ][seq_len(min(n, nrow(t))), ]
  t$symbol <- ifelse(is.na(sym[t$entrez]), t$entrez, sym[t$entrez])
  t
}

volcano <- function(tab, sym, title, subtitle, path) {
  d <- tab[!is.na(tab$padj) & !is.na(tab$lfc_shrunk), ]
  d$call <- ifelse(d$padj >= ALPHA, "not significant",
             ifelse(d$lfc_shrunk > LFC_MIN, "up (padj<0.05, LFC>1)",
             ifelse(d$lfc_shrunk < -LFC_MIN, "down (padj<0.05, LFC<-1)",
                    "significant, |LFC|<1")))
  d$call <- factor(d$call, levels = c("up (padj<0.05, LFC>1)", "down (padj<0.05, LFC<-1)",
                                      "significant, |LFC|<1", "not significant"))
  lab <- label_genes(tab, sym)

  p <- ggplot(d, aes(lfc_shrunk, -log10(pmax(padj, 1e-300)))) +
    geom_point(aes(colour = call), size = 0.5, alpha = 0.6) +
    geom_vline(xintercept = c(-LFC_MIN, LFC_MIN), linetype = 3, colour = GREY, linewidth = 0.3) +
    geom_hline(yintercept = -log10(ALPHA), linetype = 3, colour = GREY, linewidth = 0.3) +
    geom_text_repel(data = lab, aes(label = symbol), size = 2.3, colour = "black",
                    max.overlaps = 20, box.padding = 0.3, seed = 1,
                    segment.colour = GREY, segment.size = 0.25) +
    scale_colour_manual(values = setNames(c(RED, BLUE, GREY, PALE), levels(d$call)), name = NULL) +
    labs(x = "log2 fold change (apeglm shrunken)",
         y = expression(-log[10]~adjusted~italic(p)),
         title = title, subtitle = subtitle) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(hjust = 0, face = "plain"),
          plot.subtitle = element_text(hjust = 0, colour = GREY, size = 7),
          legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6.5)) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6, alpha = 1)))
  ggsave(path, plot = p, width = 5.2, height = 4.6, dpi = 300)
  path
}

ma_plot <- function(tab, sym, title, path) {
  d <- tab[!is.na(tab$padj) & !is.na(tab$lfc_shrunk), ]
  d$sig <- d$padj < ALPHA
  lab <- label_genes(tab, sym, n = 8)
  p <- ggplot(d, aes(baseMean, lfc_shrunk)) +
    geom_point(aes(colour = sig), size = 0.4, alpha = 0.5) +
    geom_hline(yintercept = 0, colour = GREY, linewidth = 0.3) +
    geom_text_repel(data = lab, aes(label = ifelse(is.na(sym[entrez]), entrez, sym[entrez])),
                    size = 2.2, colour = "black", seed = 1, max.overlaps = 15,
                    segment.colour = GREY, segment.size = 0.25) +
    scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_colour_manual(values = c(`TRUE` = RED, `FALSE` = PALE),
                        labels = c(`TRUE` = "padj < 0.05", `FALSE` = "not significant"),
                        name = NULL) +
    labs(x = "mean normalised count", y = "log2 fold change (shrunken)", title = title) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(hjust = 0, face = "plain"),
          legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6.5))
  ggsave(path, plot = p, width = 5.2, height = 3.8, dpi = 300)
  path
}
