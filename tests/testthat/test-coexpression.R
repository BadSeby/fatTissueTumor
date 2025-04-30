test_that("gene_coexpression calculates the correct correlation", {
  data("expr_small", package = "fatTissueTumor")
  norm <- normalize_expression(expr_small, method = "zscore")
  gene1 <- rownames(norm)[1]
  gene2 <- rownames(norm)[2]
  cor_val <- gene_coexpression(norm, gene1, gene2)
  expect_true(is.numeric(cor_val))
  expect_true(cor_val <= 1 && cor_val >= -1)
})
