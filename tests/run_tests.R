#!/usr/bin/env Rscript
# Run from the project root: Rscript tests/run_tests.R
PROJECT_ROOT <- normalizePath(".")
Sys.setenv(PROJECT_ROOT = PROJECT_ROOT)
library(testthat)
res <- test_file("tests/test_pipeline.R", reporter = "summary")
df <- as.data.frame(res)
failed <- sum(df$failed) + sum(df$error)
cat(sprintf("\n%d passed, %d failed/errored, %d skipped\n",
            sum(df$passed), failed, sum(df$skipped)))
quit(status = if (failed > 0) 1 else 0)
