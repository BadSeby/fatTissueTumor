#' Display a KEGG pathway with expression data
#'
#' Use pathview to display a KEGG pathway with expression data overlay.
#'
#' @param gene_data Named vector: expression values (or fold change), names = EntrezID or symbols.
#' @param pathway_id KEGG ID of pathway (e.g. "hsa04910" for Insulin signaling).
#' @param species KEGG code of the species (default: "hsa" for Homo sapiens).
#' @param out_dir Output folder for PNG files (default: tempdir()).
#' @param ... Other parameters passed to pathview.
#' @return Path to the generated PNG file.
#' @export
#' @importFrom utils data head
plot_pathway_kegg <- function(gene_data, pathway_id, species = "hsa", out_dir = tempdir(), ...) {
  if (!requireNamespace("pathview", quietly = TRUE)) {
    stop("The pathview package is required for this function. Please install it.")
  }
  # Carica bods se serve (workaround per bug pathview >= 1.41)
  if (!exists("bods", where = asNamespace("pathview"), inherits = FALSE)) {
    data("bods", package = "pathview", envir = asNamespace("pathview"))
  }
  if (is.null(names(gene_data))) stop("gene_data must have names (EntrezID or symbols).")
  res <- pathview::pathview(
    gene.data  = gene_data,
    pathway.id = pathway_id,
    species    = species,
    out.suffix = "fatTissueTumor",
    kegg.dir   = out_dir,
    ...
  )
  png_file <- list.files(out_dir, pattern = paste0(pathway_id, ".*\\.png$"), full.names = TRUE)
  if (length(png_file) == 0) warning("No PNG files found.")
  return(png_file)
}
