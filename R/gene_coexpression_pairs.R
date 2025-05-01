#' Identifies the most co-expressed gene pairs
#'
#' Calculates the correlation between all gene pairs and returns the most co-expressed.
#'
#' @param expr_data Expression matrix (genes x samples).
#' @param top_n Number of pairs to return (default: 10).
#' @param method Correlation method: 'pearson’, ‘spearman’, ‘kendall’ (default: ‘pearson’).
#' @return Data frame with most co-expressed gene pairs and correlation value.
#' @examples
#' data(expr_breast)
#' top_pairs <- top_gene_coexpression(expr_breast, top_n = 5)
#' print(top_pairs)
#' @export
top_gene_coexpression <- function(expr_data, method = "pearson", top_n = 10) {
  gene_names <- rownames(expr_data)
  cor_mat <- cor(t(expr_data), method = method)
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  res <- which(!is.na(cor_mat), arr.ind = TRUE)
  df <- data.frame(
    gene1 = gene_names[res[, 1]],
    gene2 = gene_names[res[, 2]],
    correlation = cor_mat[!is.na(cor_mat)]
  )
  df_sorted_pos <- head(df[order(-df$correlation), ], top_n)
  df_sorted_neg <- head(df[order(df$correlation), ], top_n)
  result <- rbind(df_sorted_pos, df_sorted_neg)
  rownames(result) <- NULL
  return(result)
}
