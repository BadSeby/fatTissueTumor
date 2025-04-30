test_that("top_gene_coexpression returns a data frame with the correct columns", {
  data("expr_small", package = "fatTissueTumor")
  norm <- normalize_expression(expr_small, method = "zscore")
  filt <- filter_genes(norm, min_mean = 0, min_var = 0.5)
  res <- top_gene_coexpression(filt, method = "pearson", top_n = 3)
  expect_true(is.data.frame(res))
  expect_true(all(c("gene1", "gene2", "correlation") %in% colnames(res)))
  expect_true(nrow(res) == 6) # 3 positive + 3 negative
})
