#' View and save the most co-expressed gene pairs
#'
#' Create scatterplots for each of the first N co-expressed gene pairs, with customization options.
#'
#' @param expr_data Expression matrix (genes x samples).
#' @param coexp_df Data frame produced by top_gene_coexpression().
#' @param n_plot Number of pairs to display (default: 3).
#' @param point_col Color of points (default: "darkblue").
#' @param point_cex Size of points (default: 1.5).
#' @param save_plot If TRUE, save graphs as PNG (default: FALSE).
#' @param prefix Prefix for PNG files (default: "gene_pair_").
#' @examples
#' data(expr_small)
#' filt_expr <- filter_genes(expr_small, min_mean = 1, min_var = 0.1)
#' coexp <- top_gene_coexpression(filt_expr, top_n = 5)
#' plot_top_gene_pairs(filt_expr, coexp, n_plot = 3, save_plot = TRUE)
#' @return No return value, called for side effects (plot).
#' @export
plot_top_gene_pairs <- function(expr_data, coexp_df, n_plot = 3,
                                point_col = "darkblue", point_cex = 1.5,
                                save_plot = FALSE, prefix = "gene_pair_") {
  n_plot <- min(n_plot, nrow(coexp_df))
  op <- par(mfrow = c(1, n_plot))
  for (i in seq_len(n_plot)) {
    gene1 <- coexp_df$gene1[i]
    gene2 <- coexp_df$gene2[i]
    corr <- round(coexp_df$correlation[i], 2)
    main_title <- paste0(gene1, " vs ", gene2, "\ncor=", corr)
    if (save_plot) {
      png(filename = paste0(prefix, i, ".png"), width = 700, height = 500)
    }
    plot(expr_data[gene1, ], expr_data[gene2, ],
         xlab = gene1, ylab = gene2,
         main = main_title,
         pch = 19, col = point_col, cex = point_cex)
    abline(lm(expr_data[gene2, ] ~ expr_data[gene1, ]), col = "red", lwd = 2)
    if (save_plot) dev.off()
  }
  par(op)
}
