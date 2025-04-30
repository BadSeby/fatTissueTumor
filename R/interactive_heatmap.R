#' Interactive heatmap of a gene expression
#'
#' Create an interactive heatmap of expression data using heatmaply.
#'
#' @param expr_matrix Expression matrix (genes x samples).
#' @param top_genes Number of genes plus variables to display (default: 50).
#' @param scale If TRUE, scale rows (default: TRUE).
#' @param ... Other parameters passed to heatmaply.
#' @return htmlwidget object (viewable in RStudio or browser).
#' @export
interactive_heatmap <- function(expr_matrix, top_genes = 50, scale = TRUE, ...) {
  if (!requireNamespace("heatmaply", quietly = TRUE)) {
    stop("The heatmaply package is required for this function. Please install it.")
  }
  # Seleziona i top_genes più variabili
  vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
  top_idx <- order(vars, decreasing = TRUE)[1:min(top_genes, length(vars))]
  mat <- expr_matrix[top_idx, , drop = FALSE]
  if (scale) {
    mat <- t(scale(t(mat)))
  }
  heatmaply::heatmaply(
    mat,
    scale = "none",
    xlab = "Campioni",
    ylab = "Geni",
    main = paste("Interactive heatmap of the", nrow(mat), "more variables genes"),
    ...
  )
}
