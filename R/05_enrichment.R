#!/usr/bin/env Rscript
# 05_enrichment.R — fgsea over Reactome sets, ranked by the Wald statistic.
#
# Ranking uses stat (the unshrunken Wald statistic), not the shrunken LFC:
# GSEA wants a ranking that reflects evidence strength, and stat already
# combines effect and precision. Ranking by shrunken LFC would double-count
# the shrinkage prior.
#
# Caveat that must travel with these results: every mouse Reactome annotation
# in NCBI2Reactome.txt carries the IEA evidence code (inferred electronically,
# mostly by orthology from human). None is a traceable experimental statement
# in mouse. So an enriched pathway here is a hypothesis about mouse biology
# inherited from human curation, not independent mouse evidence.

suppressPackageStartupMessages({
  library(fgsea); library(dplyr); library(ggplot2)
})

FIG <- "results/figures"
NPERM_EPS <- 0     # fgsea multilevel; eps=0 for exact small p-values
GREY <- "#5A5A5A"; BLUE <- "#1F5FA8"; RED <- "#8C2D2D"

rank_genes <- function(tab) {
  r <- tab[!is.na(tab$stat), ]
  v <- setNames(r$stat, r$entrez)
  sort(v, decreasing = TRUE)
}

run_fgsea <- function(tab, sets, min_size = 10, max_size = 500, seed = 1) {
  set.seed(seed)
  r <- rank_genes(tab)
  # fgsea prints a progress bar to stderr; it floods non-interactive logs
  res <- suppressMessages(withCallingHandlers(
    fgsea(pathways = sets, stats = r, minSize = min_size,
          maxSize = max_size, eps = NPERM_EPS, nproc = 1),
    message = function(m) invokeRestart("muffleMessage")))
  res <- res[order(res$padj), ]
  res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), "")
  as.data.frame(res)
}

plot_enrichment <- function(res, title, path, n_each = 8) {
  # Take the top n_each in EACH direction rather than the top 2*n_each by
  # |NES|. In this dataset the down side is far stronger, so a pure |NES|
  # ranking returns only negative bars and hides the up-regulated half of
  # the result entirely.
  d <- res[!is.na(res$padj) & res$padj < 0.05, ]
  up <- head(d[order(-d$NES), ], n_each)
  dn <- head(d[order( d$NES), ], n_each)
  d <- rbind(up, dn)
  d <- d[!duplicated(d$pathway), ]
  d$pathway <- factor(d$pathway, levels = d$pathway[order(d$NES)])
  d$dir <- ifelse(d$NES > 0, "enriched up", "enriched down")
  n_up <- sum(d$NES > 0); n_dn <- sum(d$NES < 0)
  wrap <- function(x) {
    # Cap at two lines: Reactome names run to 70+ characters and a free wrap
    # produces 3-4 line labels that collide with their neighbours.
    vapply(x, function(s) {
      w <- strwrap(s, width = 34)
      if (length(w) > 2) w <- c(w[1], paste0(substr(w[2], 1, 31), "..."))
      paste(w, collapse = "\n")
    }, "")
  }
  p <- ggplot(d, aes(NES, pathway, fill = dir)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
    scale_fill_manual(values = c(`enriched up` = RED, `enriched down` = BLUE), name = NULL) +
    scale_y_discrete(labels = wrap) +
    labs(x = "normalised enrichment score", y = NULL, title = title,
         subtitle = sprintf("Reactome sets at padj<0.05: top %d up and top %d down by NES (%d and %d shown).\nMouse annotations are all IEA (inferred by orthology), so these are hypotheses\ncarried over from human curation, not independent mouse evidence.",
                            n_each, n_each, n_up, n_dn)) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(hjust = 0, face = "plain"),
          plot.subtitle = element_text(hjust = 0, colour = GREY, size = 6.5),
          axis.text.y = element_text(size = 6, lineheight = 0.9),
          legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6.5))
  ggsave(path, plot = p, width = 6.6, height = 5.2, dpi = 300)
  path
}
