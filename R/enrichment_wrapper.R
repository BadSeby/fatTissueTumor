#' Universal enrichment wrapper: automatic ID conversion
#'
#' Accepts a list of genes (probeID, symbols, Entrez, Ensembl) and automatically converts to HGNC symbols.
#' Then runs enrichment_analysis().
#'
#' @param gene_list Vector of gene IDs (probe, symbols, Entrez, Ensembl).
#' @param id_type Type of ID if known ("auto", "symbol", "entrez", "ensembl", "probe"). Default: 'auto.
#' @param platform Platform name (e.g., "hgu133a", "hgu133plus2") if probe. Default: NULL.
#' @param ... Other parameters passed to enrichment_analysis().
#' @return object enrichResult.
#' @export
enrichment_auto <- function(gene_list, id_type = "auto", platform = NULL, ...) {
  # 1. Se id_type = auto, prova a indovinare
  if (id_type == "auto") {
    if (all(grepl("^ENSG", gene_list))) {
      id_type <- "ensembl"
    } else if (all(grepl("^[0-9]+$", gene_list))) {
      id_type <- "entrez"
    } else if (all(grepl("_at$", gene_list))) {
      id_type <- "probe"
    } else if (all(nchar(gene_list) <= 12)) {
      id_type <- "symbol"
    } else {
      id_type <- "symbol" # fallback
    }
  }
  # 2. Conversione in simboli genici
  if (id_type == "symbol") {
    gene_symbols <- gene_list
  } else if (id_type == "entrez") {
    # Entrez -> Symbol
    mart <- biomaRt::useEnsembl("genes", "hsapiens_gene_ensembl")
    res <- biomaRt::getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters = "entrezgene_id",
      values = gene_list,
      mart = mart
    )
    gene_symbols <- unique(na.omit(res$hgnc_symbol))
  } else if (id_type == "ensembl") {
    # Ensembl -> Symbol
    mart <- biomaRt::useEnsembl("genes", "hsapiens_gene_ensembl")
    res <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = gene_list,
      mart = mart
    )
    gene_symbols <- unique(na.omit(res$hgnc_symbol))
  } else if (id_type == "probe") {
    # Probe -> Symbol (serve platform)
    if (is.null(platform)) stop("Specify the name of the platform (e.g. hgu133a).")
    pkg <- paste0(platform, ".db")
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE)
    }
    db <- getExportedValue(pkg, pkg)
    gene_symbols <- AnnotationDbi::mapIds(
      db,
      keys = gene_list,
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )
    gene_symbols <- na.omit(gene_symbols)
  } else {
    stop("Type of ID not recognized.")
  }
  # 3. Arricchimento
  enrichment_analysis(gene_symbols, ...)
}
