#' Normalize a gene expression matrix
#'
#' Apply different normalization methods to an expression matrix (genes x samples).
#'
#' @param expr_data Expression matrix (genes x samples).
#' @param method Normalization method: "log2", "zscore", "quantile". Default: 'log2.
#' @return Normalized matrix.
#' @examples
#' data(expr_small)
#' norm_expr <- normalize_expression(expr_small, method = "zscore")
#' @export
normalize_expression <- function(expr_data, method = "log2") {
  if (method == "log2") {
    expr_norm <- log2(expr_data + 1)
  } else if (method == "zscore") {
    expr_norm <- t(scale(t(expr_data)))
  } else if (method == "quantile") {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("nstall the preprocessCore package for quantile normalization.")
    }
    expr_norm <- preprocessCore::normalize.quantiles(as.matrix(expr_data))
    rownames(expr_norm) <- rownames(expr_data)
    colnames(expr_norm) <- colnames(expr_data)
  } else {
    stop("Unrecognized method. Use 'log2', 'zscore' or 'quantile'.")
  }
  return(expr_norm)
}

#' Filters genes with low expression or low variance
#'
#' Removes genes with mean expression or variance below a threshold.
#'
#' @param expr_data Expression matrix (genes x samples).
#' @param min_mean Minimum average expression threshold (default: 1).
#' @param min_var Minimum variance threshold (default: 0.1).
#' @return Filtered matrix.
#' @examples
#' data(expr_small)
#' filt_expr <- filter_genes(expr_small, min_mean = 1, min_var = 0.1)
#' @export
filter_genes <- function(expr_data, min_mean = 1, min_var = 0.1) {
  gene_means <- rowMeans(expr_data)
  gene_vars <- apply(expr_data, 1, var)
  keep <- (gene_means >= min_mean) & (gene_vars >= min_var)
  expr_data[keep, , drop = FALSE]
}

#' Correct batch effect with ComBat
#'
#' Apply batch effect correction using ComBat.
#'
#' @param expr_data Expression matrix (genes x samples).
#' @param batch Vector indicating the batch of each sample.
#' @return Corrected matrix.
#' @examples
#' data(expr_small)
#' batch <- c(rep(1, 5), rep(2, 5))
#' expr_combat <- correct_batch(expr_small, batch)
#' @export
correct_batch <- function(expr_data, batch) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Install the sva package is required for this function. Please install it.")
  }
  expr_corrected <- sva::ComBat(dat = as.matrix(expr_data), batch = batch, par.prior = TRUE, prior.plots = FALSE)
  return(expr_corrected)
}
