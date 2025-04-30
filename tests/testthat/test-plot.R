test_that("plot_top_gene_pairs gives no errors", {
  data("expr_small", package = "fatTissueTumor")
  norm <- normalize_expression(expr_small, method = "zscore")
  filt <- filter_genes(norm, min_mean = 0, min_var = 0.5)
  res <- top_gene_coexpression(filt, method = "pearson", top_n = 3)
  expect_silent(plot_top_gene_pairs(filt, res, n_plot = 2))
})
