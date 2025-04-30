test_that("normalize_expression restituisce una matrice della stessa dimensione", {
  data("expr_small", package = "fatTissueTumor")
  norm <- normalize_expression(expr_small, method = "zscore")
  expect_equal(dim(norm), dim(expr_small))
})

test_that("filter_genes filters correctly", {
  data("expr_small", package = "fatTissueTumor")
  norm <- normalize_expression(expr_small, method = "zscore")
  filt <- filter_genes(norm, min_mean = 0, min_var = 0.5)
  expect_true(nrow(filt) <= nrow(norm))
  expect_equal(ncol(filt), ncol(norm))
})

test_that("correct_batch returns an array of the same size", {
  data("expr_small", package = "fatTissueTumor")
  norm <- normalize_expression(expr_small, method = "zscore")
  filt <- filter_genes(norm, min_mean = 0, min_var = 0.5)
  batch <- c(rep(1, 5), rep(2, 5))
  if (requireNamespace("sva", quietly = TRUE)) {
    batch_corr <- correct_batch(filt, batch)
    expect_equal(dim(batch_corr), dim(filt))
  }
})
