#' Download and import expression data from GEO (with automatic filters)
#'
#' Download a GEO series and return a list with expression matrix, phenotype table and feature data.
#'
#' @param gse_id ID of GEO series (e.g., "GSE2508").
#' @param platform_id platform ID (optional, e.g. "GPL570").
#' @param sample_type Filter on sample type (optional, e.g. "tumor", "normal").
#' @param gene_list List of genes to keep (optional).
#' @param max_na_frac Maximum fraction of NA allowed per gene (default: 0.2).
#' @param normalize If TRUE, normalize with z-score (default: FALSE).
#' @param convert_id If "symbol" or "entrez", converts probe ID to symbol or EntrezID (default: NULL).
#' @return A list with expression, phenoData, featureData.
#' @examples
#' gse <- import_GEO("GSE2508", sample_type = "obese", normalize = TRUE)
#' str(gse)
#' @export
import_GEO <- function(gse_id,
                       platform_id = NULL,
                       sample_type = NULL,
                       gene_list = NULL,
                       max_na_frac = 0.2,
                       normalize = FALSE,
                       convert_id = NULL) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("The GEOquery package is required for this function. Please install it.")
  }
  gse <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  # Se più piattaforme, scegli quella desiderata
  if (!is.null(platform_id)) {
    gse <- gse[which(sapply(gse, Biobase::annotation) == platform_id)]
    if (length(gse) == 0) stop("Platform not found.")
    gse <- gse[[1]]
  } else if (length(gse) > 1) {
    gse <- gse[[1]]
  }
  expr <- Biobase::exprs(gse)
  pheno <- Biobase::pData(gse)
  feature <- Biobase::fData(gse)
  # Filtro su tipo campione (se specificato)
  if (!is.null(sample_type)) {
    keep <- grepl(sample_type, apply(pheno, 1, paste, collapse = " "), ignore.case = TRUE)
    expr <- expr[, keep, drop = FALSE]
    pheno <- pheno[keep, , drop = FALSE]
  }
  # Filtro su NA
  na_frac <- rowMeans(is.na(expr))
  expr <- expr[na_frac <= max_na_frac, , drop = FALSE]
  feature <- feature[rownames(expr), , drop = FALSE]
  # Filtro su lista di geni
  if (!is.null(gene_list)) {
    if ("Gene Symbol" %in% colnames(feature)) {
      keep <- feature$`Gene Symbol` %in% gene_list
      expr <- expr[keep, , drop = FALSE]
      feature <- feature[keep, , drop = FALSE]
    }
  }
  # Conversione ID
  if (!is.null(convert_id)) {
    if ("Gene Symbol" %in% colnames(feature) && convert_id == "symbol") {
      rownames(expr) <- feature$`Gene Symbol`
    }
    if ("ENTREZ_GENE_ID" %in% colnames(feature) && convert_id == "entrez") {
      rownames(expr) <- feature$ENTREZ_GENE_ID
    }
  }
  # Normalizzazione opzionale
  if (normalize) {
    expr <- t(scale(t(expr)))
  }
  return(list(expr = expr, pheno = pheno, feature = feature))
}

#' Download and import expression data from TCGA (with automatic filters)
#'
#' Download RNA-seq data from TCGA and return a list with expression matrix and phenotype table.
#'
#' @param project TCGA project name (e.g., "TCGA-BRCA").
#' @param sample_type Filter on sample type (e.g., "Primary Tumor," "Solid Tissue Normal").
#' @param gene_list List of genes to keep (optional).
#' @param normalize If TRUE, normalize with z-score (default: FALSE).
#' @return A list with expression and phenoData.
#' @examples
#' \dontrun{
#' tcga <- import_TCGA("TCGA-BRCA", sample_type = "Primary Tumor",
#'  normalize = TRUE, workflow.type = "STAR - Counts")
#' str(tcga)
#' }
#' @export
import_TCGA <- function(project,
                        sample_type = NULL,
                        gene_list = NULL,
                        normalize = FALSE) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("The TCGAbiolinks package is required for this function. Please install it.")
  }
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts"
  )
  TCGAbiolinks::GDCdownload(query)
  data <- TCGAbiolinks::GDCprepare(query)
  expr <- SummarizedExperiment::assay(data)
  pheno <- as.data.frame(SummarizedExperiment::colData(data))
  # Filtro su tipo campione
  if (!is.null(sample_type)) {
    keep <- pheno$short_letter_code %in% sample_type | pheno$sample_type %in% sample_type
    expr <- expr[, keep, drop = FALSE]
    pheno <- pheno[keep, , drop = FALSE]
  }
  # Filtro su lista di geni
  if (!is.null(gene_list)) {
    expr <- expr[rownames(expr) %in% gene_list, , drop = FALSE]
  }
  # Normalizzazione opzionale
  if (normalize) {
    expr <- t(scale(t(expr)))
  }
  return(list(expr = expr, pheno = pheno))
}

#' Download and import expression data from ArrayExpress (with automatic filters)
#'
#' Download an ArrayExpress experiment and return a list with expression matrix, phenotype table and feature data.
#'
#' @param ae_id ArrayExpress experiment ID (e.g., "E-MTAB-62").
#' @param gene_list List of genes to keep (optional).
#' @param max_na_frac Maximum fraction of NA allowed per gene (default: 0.2).
#' @param normalize If TRUE, normalize with z-score (default: FALSE).
#' @return A list with expression, phenoData, featureData.
#' @examples
#' \dontrun{
#' # This example requires an internet connection and may fail during automated checks
#' ae <- import_ArrayExpress("E-MTAB-62", normalize = TRUE)
#' str(ae)
#' }
#' @export
import_ArrayExpress <- function(ae_id,
                                gene_list = NULL,
                                max_na_frac = 0.2,
                                normalize = FALSE) {
  if (!requireNamespace("ArrayExpress", quietly = TRUE)) {
    stop("The ArrayExpress package is required for this function. Please install it.")
  }
  ae <- ArrayExpress::ArrayExpress(ae_id)
  expr <- Biobase::exprs(ae)
  pheno <- Biobase::pData(ae)
  feature <- Biobase::fData(ae)
  # Filtro su NA
  na_frac <- rowMeans(is.na(expr))
  expr <- expr[na_frac <= max_na_frac, , drop = FALSE]
  feature <- feature[rownames(expr), , drop = FALSE]
  # Filtro su lista di geni
  if (!is.null(gene_list)) {
    if ("Gene Symbol" %in% colnames(feature)) {
      keep <- feature$`Gene Symbol` %in% gene_list
      expr <- expr[keep, , drop = FALSE]
      feature <- feature[keep, , drop = FALSE]
    }
  }
  # Normalizzazione opzionale
  if (normalize) {
    expr <- t(scale(t(expr)))
  }
  return(list(expr = expr, pheno = pheno, feature = feature))
}
