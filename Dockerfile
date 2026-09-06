# Pinned Bioconductor base: R and every Bioconductor package version is fixed
# by the tag, so the image is the reproducibility contract, not the R code.
FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL org.opencontainers.image.title="GSE60450 RNA-seq differential expression"
LABEL org.opencontainers.image.description="DESeq2 + fgsea pipeline with a composition-bias sensitivity analysis"

WORKDIR /pipeline

# CRAN packages first (slowest layer, changes least)
RUN R -q -e "install.packages(c('tidyverse','ggrepel','pheatmap','rmarkdown','knitr','testthat','sessioninfo','jsonlite'), repos='https://cloud.r-project.org')"

# Bioconductor packages. org.Mm.eg.db is deliberately NOT installed: the
# pipeline reads gene symbols from NCBI gene_info and pathways from Reactome,
# both fetched by download_data.sh, so annotation provenance is a file rather
# than a package snapshot.
RUN R -q -e "BiocManager::install(c('DESeq2','apeglm','fgsea'), ask=FALSE, update=FALSE)"

COPY R/ ./R/
COPY tests/ ./tests/
COPY run_pipeline.R build_sample_sheet.R download_data.sh ./
RUN chmod +x download_data.sh

# Inputs are not baked in — they are public and ~40 MB. Run:
#   docker run --rm -v "$PWD/results:/pipeline/results" IMAGE \
#     bash -c "./download_data.sh && Rscript run_pipeline.R"
CMD ["bash", "-c", "./download_data.sh && Rscript tests/run_tests.R && Rscript run_pipeline.R"]
