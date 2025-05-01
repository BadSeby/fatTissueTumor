#' GO, KEGG, Reactome enrichment analysis
#'
#' Performs functional enrichment on a list of genes.
#'
#' @param gene_list Gene symbol vector or EntrezID.
#' @param organism Organism (default: "hsa" for Homo sapiens).
#' @param pvalueCutoff Significance threshold (default: 0.05).
#' @param type Type of analysis: GO, KEGG, Reactome.
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   # Example gene symbols (human)
#'   gene_list <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN")
#'   # GO enrichment analysis
#'   ego <- enrichment_analysis(gene_list, type = "GO")
#'   if (methods::is(ego, "enrichResult") && nrow(ego@result) > 0) {
#'     print(head(ego@result))
#'   }
#'   # KEGG enrichment analysis (requires Entrez IDs)
#'   # Convert gene symbols to Entrez IDs
#'   if (requireNamespace("AnnotationDbi", quietly = TRUE)) {
#'     entrez_ids <- AnnotationDbi::mapIds(
#'       org.Hs.eg.db::org.Hs.eg.db,
#'       keys = gene_list,
#'       column = "ENTREZID",
#'       keytype = "SYMBOL",
#'       multiVals = "first"
#'     )
#'     entrez_ids <- na.omit(entrez_ids)
#'     if (length(entrez_ids) > 0) {
#'       ekegg <- enrichment_analysis(unname(entrez_ids), type = "KEGG")
#'       if (methods::is(ekegg, "enrichResult") && nrow(ekegg@result) > 0) {
#'         print(head(ekegg@result))
#'       }
#'     }
#'   }
#'   # Reactome enrichment analysis (if ReactomePA is available)
#'   if (requireNamespace("ReactomePA", quietly = TRUE)) {
#'     ereactome <- enrichment_analysis(gene_list, type = "Reactome")
#'     if (methods::is(ereactome, "enrichResult") && nrow(ereactome@result) > 0) {
#'       print(head(ereactome@result))
#'     }
#'   }
#' }
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
#' @examples
#' \dontrun{
#' if (requireNamespace("enrichR", quietly = TRUE)) {
#'   # Example with genes symbols
#'   gene_list <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN")
#'   dbs <- c("KEGG_2021_Human")
#'   results <- enrichment_enrichr(gene_list, databases = dbs)
#'   # View first rows
#'   if (length(results) > 0) print(head(results[[1]]))
#' }
#' }
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
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("enrichplot", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   gene <- c("4312", "8318", "10874", "55143", "55388", "991", "6280", "23046", "949", "10533")
#'   ego <- clusterProfiler::enrichGO(
#'     gene          = gene,
#'     OrgDb         = "org.Hs.eg.db",
#'     keyType       = "ENTREZID",
#'     ont           = "BP",
#'     pAdjustMethod = "BH",
#'     pvalueCutoff  = 0.1,
#'     readable      = FALSE
#'   )
#'   if (nrow(ego@result) > 0) {
#'     plot_enrichment(ego, top = 5)
#'   }
#' }
#' @export
plot_enrichment <- function(enrich_result, top = 10) {
  barplot(enrich_result, showCategory = top)
}


#' Dotplot of enrichment results
#'
#' @param enrich_result Object enrichResult.
#' @param top No. of terms to display (default: 10).
#' @return No return value, called for side effects (plot).
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("enrichplot", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   # Example with Entrez IDs (no internet required)
#'   gene <- c("4312", "8318", "10874", "55143", "55388", "991", "6280", "23046", "949", "10533")
#'   ego <- clusterProfiler::enrichGO(
#'     gene          = gene,
#'     OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
#'     keyType       = "ENTREZID",
#'     ont           = "BP",
#'     pAdjustMethod = "BH",
#'     pvalueCutoff  = 0.1,
#'     readable      = FALSE
#'   )
#'   if (methods::is(ego, "enrichResult") && nrow(ego@result) > 0) {
#'     dotplot_enrichment(ego, top = 5)
#'   }
#' }
#' @export
dotplot_enrichment <- function(enrich_result, top = 10) {
  dotplot(enrich_result, showCategory = top)
}
