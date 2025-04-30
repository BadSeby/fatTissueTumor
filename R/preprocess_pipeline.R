#' Advanced preprocessing pipeline for expression data
#'
#' Performs sequentially: missing removal, filters on genes and samples, normalization, log-transform, batch correction.
#'
#' @param expr_matrix Expression matrix (genes x samples).
#' @param min_expr Minimum average expression value to hold a gene (default: 1).
#' @param min_expr_frac Minimum fraction of samples in which gene must exceed min_expr (default: 0.2).
#' @param min_var Minimum variance threshold to hold a gene (default: 0.01).
#' @param max_na_frac Maximum fraction of NA allowed per gene (default: 0.2).
#' @param min_expr_sample Minimum mean expression value to hold a sample (default: 1).
#' @param min_var_sample Minimum variance threshold to hold a sample (default: 0).
#' @param max_na_frac_sample Maximum fraction of NA allowed per sample (default: 0.2).
#' @param normalize Normalization method: 'none’, ‘zscore’, ‘quantile’, ‘minmax’ (default: ‘zscore’).
#' @param log_transform If TRUE, apply log2(x+1) (default: FALSE).
#' @param batch Batch vector for batch correction (optional).
#' @param return_report If TRUE, also returns a detailed report (default: FALSE).
#' @return Preprocessed matrix (and, if required, a report).
#' @examples
#' \dontrun{
#' data(expr_breast)
#' res <- preprocess_expression(expr_breast, log_transform = TRUE, return_report = TRUE)
#' str(res)
#' }
#' @export
preprocess_expression <- function(expr_matrix,
                                  min_expr = 1,
                                  min_expr_frac = 0.2,
                                  min_var = 0.01,
                                  max_na_frac = 0.2,
                                  min_expr_sample = 1,
                                  min_var_sample = 0,
                                  max_na_frac_sample = 0.2,
                                  normalize = "zscore",
                                  log_transform = FALSE,
                                  batch = NULL,
                                  return_report = FALSE) {
  mat <- as.matrix(expr_matrix)
  report <- list()
  report$start_dim <- dim(mat)

  # 1. Rimozione geni con troppi NA
  na_frac <- rowMeans(is.na(mat))
  keep_gene_na <- na_frac <= max_na_frac
  report$genes_na_removed <- sum(!keep_gene_na)
  mat <- mat[keep_gene_na, , drop = FALSE]

  # 2. Filtro per espressione minima in almeno min_expr_frac dei campioni
  keep_gene_expr <- apply(mat, 1, function(x) mean(x > min_expr, na.rm = TRUE) >= min_expr_frac)
  report$genes_expr_removed <- sum(!keep_gene_expr)
  mat <- mat[keep_gene_expr, , drop = FALSE]

  # 3. Filtro per varianza minima sui geni
  keep_gene_var <- apply(mat, 1, function(x) var(x, na.rm = TRUE) >= min_var)
  report$genes_var_removed <- sum(!keep_gene_var)
  mat <- mat[keep_gene_var, , drop = FALSE]

  # 4. Rimozione campioni con troppi NA
  na_frac_sample <- colMeans(is.na(mat))
  keep_sample_na <- na_frac_sample <= max_na_frac_sample
  report$samples_na_removed <- sum(!keep_sample_na)
  mat <- mat[, keep_sample_na, drop = FALSE]

  # 5. Filtro per espressione minima sui campioni
  keep_sample_expr <- apply(mat, 2, function(x) mean(x, na.rm = TRUE) >= min_expr_sample)
  report$samples_expr_removed <- sum(!keep_sample_expr)
  mat <- mat[, keep_sample_expr, drop = FALSE]

  # 6. Filtro per varianza minima sui campioni
  keep_sample_var <- apply(mat, 2, function(x) var(x, na.rm = TRUE) >= min_var_sample)
  report$samples_var_removed <- sum(!keep_sample_var)
  mat <- mat[, keep_sample_var, drop = FALSE]

  report$final_dim <- dim(mat)

  # 7. Log-transform
  if (log_transform) {
    mat <- log2(mat + 1)
    report$log_transform <- TRUE
  } else {
    report$log_transform <- FALSE
  }

  # 8. Normalizzazione
  if (normalize == "zscore") {
    mat <- t(scale(t(mat)))
    report$normalization <- "zscore"
  } else if (normalize == "quantile") {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("The preprocessCore package is required for quantile normalization.")
    }
    mat <- preprocessCore::normalize.quantiles(mat)
    rownames(mat) <- rownames(expr_matrix)
    colnames(mat) <- colnames(expr_matrix)
    report$normalization <- "quantile"
  } else if (normalize == "minmax") {
    mat <- t(apply(mat, 1, function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))
    report$normalization <- "minmax"
  } else {
    report$normalization <- "none"
  }

  # 9. Batch correction (ComBat)
  if (!is.null(batch)) {
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("The sva package is required for this function. Please install it.")
    }
    if (length(batch) != ncol(mat)) stop("The batch length should be equal to the number of samples.")
    mat <- sva::ComBat(dat = mat, batch = batch, par.prior = TRUE, prior.plots = FALSE)
    report$batch_correction <- TRUE
  } else {
    report$batch_correction <- FALSE
  }

  if (return_report) {
    return(list(expr = mat, report = report))
  } else {
    return(mat)
  }
}
