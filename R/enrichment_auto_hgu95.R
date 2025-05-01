#' Automatic enrichment for probe ID Affymetrix HG-U95 (A/B/C/D/E)
#'
#' Automatically map probe IDs of HG-U95A/B/C/D/E platforms to gene symbols and launch enrichment_analysis().
#'
#' @param probe_list Vector of probe IDs (e.g. "1000_at", "1600_s_at", etc.)
#' @param ... Other parameters passed to enrichment_analysis (e.g., type = "GO")
#' @return Object enrichResult
#' @export
enrichment_auto_hgu95 <- function(probe_list, ...) {
  platforms <- c("hgu95a.db", "hgu95b.db", "hgu95c.db", "hgu95d.db", "hgu95e.db")
  all_symbols <- character(0)
  for (pkg in platforms) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package HG-U95 is required for this function. Please install it.")
    }
    db <- getExportedValue(pkg, sub(".db$", ".db", pkg))
    valid_probes <- intersect(probe_list, AnnotationDbi::keys(db, keytype = "PROBEID"))
    if (length(valid_probes) > 0) {
      symbols <- AnnotationDbi::mapIds(
        db,
        keys = valid_probes,
        column = "SYMBOL",
        keytype = "PROBEID",
        multiVals = "first"
      )
      all_symbols <- c(all_symbols, na.omit(symbols))
    }
  }
  all_symbols <- unique(all_symbols)
  if (length(all_symbols) == 0) {
    stop("No probe ID mapped to any HG-U95 platform.")
  }
  enrichment_analysis(all_symbols, ...)
}
