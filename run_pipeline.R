#!/usr/bin/env Rscript
# run_pipeline.R — end-to-end: counts -> DE -> enrichment -> sensitivity.
# Usage:  Rscript run_pipeline.R
# Inputs: data/raw/*.gz, data/sample_sheet.csv, data/annotation/*
# Writes: results/*.csv, results/*.rds, results/figures/*.png, results/session_info.txt

t0 <- Sys.time()
options(warn = 1)
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

step <- function(msg) cat(sprintf("[%5.1fs] %s\n", as.numeric(Sys.time() - t0, units = "secs"), msg))

step("loading counts and sample sheet")
dds <- build_dds()
saveRDS(dds, file.path(RES, "dds_raw.rds"))

step("fitting DESeq2")
dds <- run_de(dds)
saveRDS(dds, file.path(RES, "dds_fitted.rds"))

step("annotation from NCBI + Reactome")
sym  <- load_gene_symbols()
desc <- load_gene_descriptions()
sets <- load_reactome(rownames(dds))
saveRDS(list(sym = sym, desc = desc, sets = sets), file.path(RES, "annotation.rds"))
write.csv(annotation_provenance(), file.path(RES, "annotation_provenance.csv"), row.names = FALSE)

step("QC figures")
qc <- qc_table(dds)
write.csv(qc, file.path(RES, "qc_table.csv"), row.names = FALSE)
plot_libsize_vs_sizefactor(qc)
pca <- plot_pca(dds)

step("differential expression, 2 contrasts")
CONTRASTS <- list(
  B_vs_L_virgin        = c("basal_virgin", "luminal_virgin"),
  lact_vs_preg_luminal = c("luminal_lactating", "luminal_pregnant"))
tabs <- lapply(names(CONTRASTS), function(nm)
  contrast_table(dds, CONTRASTS[[nm]][1], CONTRASTS[[nm]][2], nm))
names(tabs) <- names(CONTRASTS)
saveRDS(tabs, file.path(RES, "de_tables.rds"))

summ <- do.call(rbind, lapply(names(tabs), function(nm) {
  s <- summarise_contrast(tabs[[nm]])
  data.frame(contrast = nm, n_tested = s$n_tested, n_sig = s$n_sig,
             n_strong = s$n_strong, n_strong_unshrunk = s$n_strong_unshrunk,
             n_up = s$n_up, n_down = s$n_down)
}))
write.csv(summ, file.path(RES, "de_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

step("DE figures")
for (nm in names(tabs)) {
  readr::write_csv(tabs[[nm]], file.path(RES, paste0("de_", nm, ".csv")))
  s <- summarise_contrast(tabs[[nm]])
  sub <- sprintf("%d of %d tested genes significant at padj<0.05; %d with |LFC|>1 after shrinkage. n=2 per group.",
                 s$n_sig, s$n_tested, s$n_strong)
  ttl <- if (nm == "B_vs_L_virgin") "Basal vs luminal cells, virgin mammary gland"
         else "Lactating vs pregnant luminal cells"
  volcano(tabs[[nm]], sym, ttl, sub, file.path(FIG, paste0("volcano_", nm, ".png")))
  ma_plot(tabs[[nm]], sym, ttl, file.path(FIG, paste0("ma_", nm, ".png")))
}

step("Reactome enrichment")
gsea <- lapply(tabs, function(t) run_fgsea(t, sets))
saveRDS(gsea, file.path(RES, "gsea.rds"))
for (nm in names(gsea)) {
  readr::write_csv(gsea[[nm]], file.path(RES, paste0("gsea_", nm, ".csv")))
  ttl <- if (nm == "B_vs_L_virgin") "Basal vs luminal cells, virgin mammary gland"
         else "Lactating vs pregnant luminal cells"
  plot_enrichment(gsea[[nm]], ttl, file.path(FIG, paste0("gsea_", nm, ".png")))
}

step("sensitivity: refit without milk-protein genes")
dds_raw <- readRDS(file.path(RES, "dds_raw.rds"))
milk <- milk_gene_ids(sym, rownames(dds_raw))
dds_nm <- refit_without(dds_raw, milk)
tab_nm <- contrast_table(dds_nm, "luminal_lactating", "luminal_pregnant", "lact_vs_preg_nomilk")
g_nm <- run_fgsea(tab_nm, sets)
readr::write_csv(g_nm, file.path(RES, "gsea_lact_vs_preg_nomilk.csv"))
cmp <- compare_nes(gsea$lact_vs_preg_luminal, g_nm,
                   "ribosom|translation|rRNA|Nonsense Mediated|40S|60S")
readr::write_csv(cmp, file.path(RES, "sensitivity_nomilk.csv"))
plot_sensitivity(cmp)
cat(sprintf("  %d milk genes excluded; NES correlation %.4f; max shift %.2f\n",
            length(milk), cor(cmp$NES_full, cmp$NES_nomilk, use = "complete.obs"),
            max(abs(cmp$NES_full - cmp$NES_nomilk), na.rm = TRUE)))

step("session info")
writeLines(capture.output(sessioninfo::session_info()), file.path(RES, "session_info.txt"))
step(sprintf("done — %d figures, %d tables in %s/",
             length(list.files(FIG, pattern = "png")),
             length(list.files(RES, pattern = "csv")), RES))
