#' Display a gene pair interactively
#'
#' Use plotly for an interactive scatterplot of a gene pair.
#' @param expr_data Expression matrix (genes x samples).
#' @param gene1 First gene.
#' @param gene2 Second gene.
#' @examples
#' data(expr_small)
#' filt_expr <- filter_genes(expr_small, min_mean = 1, min_var = 0.1)
#' coexp <- top_gene_coexpression(filt_expr, top_n = 5)
#' plot_gene_pair_interactive(filt_expr, coexp$gene1[1], coexp$gene2[1])
#' @export
#' @importFrom magrittr %>%
plot_gene_pair_interactive <- function(expr_data, gene1, gene2) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Install the plotly package for interactive visualization.")
  }
  plotly::plot_ly(
    x = expr_data[gene1, ],
    y = expr_data[gene2, ],
    type = "scatter",
    mode = "markers",
    marker = list(color = "darkblue", size = 12)
  ) %>%
    plotly::layout(
      title = paste(gene1, "vs", gene2),
      xaxis = list(title = gene1),
      yaxis = list(title = gene2)
    )
}
