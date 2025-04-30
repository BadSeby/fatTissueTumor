test_that("Complete Pipeline: preprocessing, co-espressione, enrichment", {
  # Loading Data
  data("expr_small", package = "fatTissueTumor")
  expect_true(is.matrix(expr_small) || is.data.frame(expr_small))

  # Normalization
  expr_z <- normalize_expression(expr_small, method = "zscore")
  expect_equal(dim(expr_z), dim(expr_small))

  # Filtering
  expr_filt <- filter_genes(expr_z, min_mean = 0, min_var = 0.5)
  expect_true(nrow(expr_filt) < nrow(expr_z))

  #  Batch-correction (simulated)
  batch <- rep(1:2, each = ncol(expr_filt)/2)
  if (requireNamespace("sva", quietly = TRUE)) {
    expr_batch <- correct_batch(expr_filt, batch)
    expect_equal(dim(expr_batch), dim(expr_filt))
  }

  # Co-expression between two genes
  gene1 <- rownames(expr_filt)[1]
  gene2 <- rownames(expr_filt)[2]
  cor_val <- gene_coexpression(expr_filt, gene1, gene2)
  expect_true(is.numeric(cor_val))
  expect_true(cor_val <= 1 && cor_val >= -1)

  # Most co-expressed gene pairs
  coexp <- top_gene_coexpression(expr_filt, method = "pearson", top_n = 5)
  expect_true(is.data.frame(coexp))
  expect_true(nrow(coexp) == 10)

  # Automatic GO Enrichment for Affymetrix HG-U95 ID Probe
  gene_list <- rownames(expr_filt)[1:40]
  if (exists("enrichment_auto_hgu95")) {
    go_res <- enrichment_auto_hgu95(gene_list, type = "GO")
    expect_true("enrichResult" %in% class(go_res))
    expect_true(nrow(as.data.frame(go_res)) > 0)
  }

  # GO enrichment by gene symbols
  gene_symbols <- c("LEP", "STAT3", "TP53", "BRCA1")
  go_res2 <- enrichment_auto(gene_symbols, id_type = "auto", type = "GO")
  expect_true("enrichResult" %in% class(go_res2))
})
