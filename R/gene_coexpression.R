#' Calculates the co-expression correlation between two genes
#'
#' This function calculates the expression correlation between two genes in a gene expression dataset.
#'
#' @param expr_data A data frame or matrix with genes as rows and samples as columns.
#' @param gene1 Name of the first gene (string).
#' @param gene2 Name of the second gene (string).
#' @param method Correlation method: pearson, spearman or kendall. Default: pearson.
#'
#' @return A numerical value representing the correlation between the two genes.
#' @examples
#' # Example with dummy data
#' expr <- data.frame(
#' geneA = c(1, 2, 3, 4, 5),
#' geneB = c(2, 4, 6, 8, 10),
#' geneC = c(5, 4, 3, 2, 1)
#' )
#' gene_coexpression(t(expr), "geneA", "geneB")
#' @export
gene_coexpression <- function(expr_data, gene1, gene2, method = "pearson") {
  if (!(gene1 %in% rownames(expr_data))) {
    stop(paste("Gene", gene1, "not found in expr_data rows"))
  }
  if (!(gene2 %in% rownames(expr_data))) {
    stop(paste("Gene", gene2, "not found in expr_data rows"))
  }
  x <- as.numeric(expr_data[gene1, ])
  y <- as.numeric(expr_data[gene2, ])
  cor(x, y, method = method)
}
