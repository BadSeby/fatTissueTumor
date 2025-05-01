#' Creates a gene co-expression network
#'
#' Constructs a co-expression graph between the most variable genes.
#'
#' @param expr_matrix Expression matrix (genes x samples).
#' @param top_genes Number of genes plus variables to consider (default: 30).
#' @param method Correlation method ("pearson" or "spearman").
#' @param threshold Absolute correlation thr to create an arc (default: 0.7).
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   set.seed(1)
#'   expr_matrix <- matrix(rnorm(100), nrow = 10)
#'   rownames(expr_matrix) <- paste0("Gene", 1:10)
#'   g <- coexpression_network(expr_matrix, top_genes = 5, threshold = 0.2)
#' }
#' @return igraph object.
#' @export
coexpression_network <- function(expr_matrix, top_genes = 30, method = "pearson", threshold = 0.7) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The igraph package is required for this function. Please install it.")
  }
  # Seleziona i top_genes più variabili
  vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
  top_idx <- order(vars, decreasing = TRUE)[1:min(top_genes, length(vars))]
  mat <- expr_matrix[top_idx, , drop = FALSE]
  # Calcola la matrice di correlazione
  cor_mat <- cor(t(mat), method = method, use = "pairwise.complete.obs")
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  edges <- which(abs(cor_mat) >= threshold, arr.ind = TRUE)
  if (nrow(edges) == 0) stop("No pair of genes exceeds the correlation threshold.")
  edge_list <- data.frame(
    from = rownames(cor_mat)[edges[,1]],
    to = colnames(cor_mat)[edges[,2]],
    weight = cor_mat[edges]
  )
  g <- igraph::graph_from_data_frame(edge_list, directed = FALSE)
  return(g)
}


#' Display a co-expression network (static, custom)
#'
#' @param g Object igraph.
#' @param palette Color function for groups (default: rainbow).
#' @param node_anno Annotation vector/factor for nodes (optional).
#' @param ... Other parameters for plot.igraph.
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   # Expression matrix
#'   set.seed(123)
#'   expr_matrix <- matrix(rnorm(50), nrow = 5)
#'   rownames(expr_matrix) <- paste0("Gene", 1:5)
#'   colnames(expr_matrix) <- paste0("Sample", 1:10)
#'   # Nwetwok co-expression
#'   cor_mat <- cor(t(expr_matrix))
#'   cor_mat[lower.tri(cor_mat, diag = TRUE)] <- 0
#'   edges <- which(abs(cor_mat) > 0.5, arr.ind = TRUE)
#'   edge_list <- data.frame(
#'     from = rownames(cor_mat)[edges[,1]],
#'     to = colnames(cor_mat)[edges[,2]],
#'     weight = cor_mat[edges]
#'   )
#'   g <- igraph::graph_from_data_frame(edge_list, directed = FALSE)
#'   # Example without annotation
#'   plot_coexpression_network(g)
#'   # Example with nodes annotation
#'   node_anno <- factor(c("A", "A", "B", "B", "C"))
#'   names(node_anno) <- igraph::V(g)$name
#'   plot_coexpression_network(g, node_anno = node_anno)
#' }
#' @export
plot_coexpression_network <- function(g, node_anno = NULL, palette = rainbow, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The igraph package is required for this function. Please install it.")
  }
  # Degree per dimensione nodi
  deg <- igraph::degree(g)
  vsize <- 10 + 5 * (deg - min(deg)) / (max(deg) - min(deg) + 1e-6)
  # Colore nodi per gruppo
  vcol <- "skyblue"
  legend_items <- NULL
  if (!is.null(node_anno)) {
    # Associa i colori ai gruppi
    groups <- as.factor(node_anno[igraph::V(g)$name])
    pal <- palette(length(levels(groups)))
    vcol <- pal[as.integer(groups)]
    legend_items <- list(groups = levels(groups), colors = pal)
  }
  # Layout
  lay <- igraph::layout_with_fr(g, weights = abs(igraph::E(g)$weight))
  plot(
    g,
    layout = lay,
    vertex.label.cex = 0.8,
    vertex.size = vsize,
    vertex.color = vcol,
    edge.width = abs(igraph::E(g)$weight) * 5,
    edge.color = ifelse(igraph::E(g)$weight > 0, "tomato", "royalblue"),
    ...
  )
  # Legenda gruppi
  if (!is.null(legend_items)) {
    legend("topleft", legend = legend_items$groups, col = legend_items$colors, pch = 19, pt.cex = 2, bty = "n")
  }
}


#' Display a co-expression network interactively
#'
#' @param g Object igraph.
#' @param node_anno Annotation vector/factor for nodes (optional).
#' @param palette Color function for groups (default: rainbow).
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE) && requireNamespace("visNetwork", quietly = TRUE)) {
#'   set.seed(1)
#'   expr_matrix <- matrix(rnorm(100), nrow = 10)
#'   rownames(expr_matrix) <- paste0("Gene", 1:10)
#'   g <- coexpression_network(expr_matrix, top_genes = 5,threshold = 0.2)
#'   plot_coexpression_network_interactive(g)
#' }
#' @export
#' @importFrom stats cor var lm na.omit
#' @importFrom grDevices dev.off png rainbow
plot_coexpression_network_interactive <- function(g, node_anno = NULL, palette = rainbow) {
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("The visNetwork package is required for this function. Please install it.")
  }
  deg <- igraph::degree(g)
  nodes <- data.frame(
    id = igraph::V(g)$name,
    label = igraph::V(g)$name,
    value = deg
  )
  if (!is.null(node_anno)) {
    groups <- as.factor(node_anno[igraph::V(g)$name])
    pal <- palette(length(levels(groups)))
    nodes$group <- groups
    nodes$color.background <- pal[as.integer(groups)]
  }
  edges <- as.data.frame(igraph::as_data_frame(g, what = "edges"))
  edges$color <- ifelse(edges$weight > 0, "tomato", "royalblue")
  visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visEdges(smooth = FALSE) %>%
    visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
}
