#' GO, KEGG, Reactome enrichment analysis
#'
#' Performs functional enrichment on a list of genes.
#'
#' @param gene_list Gene symbol vector or EntrezID.
#' @param organism Organism (default: "hsa" for Homo sapiens).
#' @param pvalueCutoff Significance threshold (default: 0.05).
#' @param type Type of analysis: GO, KEGG, Reactome.
#' @return Object of class enrichResult or similar.
#' @export
enrichment_analysis <- function(gene_list, type = "GO", organism = "hsa", pvalueCutoff = 0.05) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Install the clusterProfiler package is required for this function. Please install it.")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Install the package org.Hs.eg.db is required for this function. Please install it.")
  }
  if (type == "GO") {
    res <- clusterProfiler::enrichGO(
      gene         = gene_list,
      OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
      keyType      = "SYMBOL",
      ont          = "BP",
      pvalueCutoff = pvalueCutoff
    )
  } else if (type == "KEGG") {
    res <- clusterProfiler::enrichKEGG(
      gene         = gene_list,
      organism     = organism,
      keyType      = "kegg",
      pvalueCutoff = pvalueCutoff
    )
  } else if (type == "Reactome") {
    if (!requireNamespace("ReactomePA", quietly = TRUE)) {
      stop("Install the ReactomePA package is required for this function. Please install it.")
    }
    res <- ReactomePA::enrichPathway(
      gene         = gene_list,
      organism     = "human",
      pvalueCutoff = pvalueCutoff,
      readable     = TRUE
    )
  } else {
    stop("Type not supported. Use GO, KEGG or Reactome.")
  }
  return(res)
}


#' Enrichment analysis with enrichR
#'
#' @param gene_list Vector of gene symbols.
#' @param databases Vector of enrichR database names (e.g., "KEGG_2021_Human").
#' @return List of data frames with results.
#' @export
enrichment_enrichr <- function(gene_list, databases = c("KEGG_2021_Human", "GO_Biological_Process_2021")) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop("Install the enrichR package is required for this function. Please install it.")
  }
  enrichR::enrichr(gene_list, databases)
}


#' Display enrichment results (barplot)
#'
#' @param enrich_result Object enrichResult (from clusterProfiler or ReactomePA).
#' @param top No. of terms to display (default: 10).
#' @export
plot_enrichment <- function(enrich_result, top = 10) {
  barplot(enrich_result, showCategory = top)
}


#' Dotplot of enrichment results
#'
#' @param enrich_result Object enrichResult.
#' @param top No. of terms to display (default: 10).
#' @export
dotplot_enrichment <- function(enrich_result, top = 10) {
  dotplot(enrich_result, showCategory = top)
}
