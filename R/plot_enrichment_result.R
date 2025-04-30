#' Display enrichment results (GO, KEGG, Reactome)
#'
#' Wrapper for main enrichplot/clusterProfiler plots.
#'
#' @param enrich_obj Object enrichResult (e.g., from enrichGO, enrichKEGG, enrichPathway).
#' @param type Type of plot: 'dotplot’, ‘barplot’, ‘cnetplot’, ‘emapplot’, ‘upsetplot’.
#' @param showCategory Number of categories to show (default: 20).
#' @param ... Other parameters for the plotting function.
#' @return A ggplot2 or enrichplot object.
#' @examples
#' \dontrun{
#' # Example (requires enrichResult object)
#' # plot_enrichment_result(ego, type = "dotplot", showCategory = 10)
#' }
#' @export
#' @importFrom graphics barplot legend par abline
#' @importFrom enrichplot dotplot
plot_enrichment_result <- function(enrich_obj, type = "dotplot", showCategory = 20, ...) {
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("The enrichplot package is required for this function. Please install it.")
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("The clusterProfiler package is required for this function. Please install it.")
  }
  switch(type,
         dotplot = enrichplot::dotplot(enrich_obj, showCategory = showCategory, ...),
         barplot = barplot(enrich_obj, showCategory = showCategory, ...),
         cnetplot = enrichplot::cnetplot(enrich_obj, showCategory = showCategory, ...),
         emapplot = enrichplot::emapplot(enrich_obj, showCategory = showCategory, ...),
         upsetplot = enrichplot::upsetplot(enrich_obj, ...),
         stop("Plot type not supported.")
  )
}
