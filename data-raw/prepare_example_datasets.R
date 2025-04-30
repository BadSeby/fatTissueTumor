# data-raw/prepare_example_datasets.R

rm(list = ls())
# Installa i pacchetti se mancano
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("ArrayExpress", quietly = TRUE)) BiocManager::install("ArrayExpress")

library(GEOquery)
library(ArrayExpress)

# Funzione per creare un subset casuale di una matrice
subset_matrix <- function(mat, n_genes = 30, n_samples = 10) {
  set.seed(123)
  genes <- sample(rownames(mat), min(n_genes, nrow(mat)))
  samples <- sample(colnames(mat), min(n_samples, ncol(mat)))
  mat[genes, samples, drop = FALSE]
}

# 1. GSE2508 (adipose/breast)
gse <- getGEO("GSE2508", GSEMatrix = TRUE)[[1]]
expr_breast <- exprs(gse)
pheno_breast <- pData(gse)
expr_breast <- subset_matrix(expr_breast, n_genes = 30, n_samples = 10)
pheno_breast <- pheno_breast[colnames(expr_breast), , drop = FALSE]
save(expr_breast, pheno_breast, file = "expr_breast.rda")

# 2. GSE44076 (colon)
gse2 <- getGEO("GSE44076", GSEMatrix = TRUE)[[1]]
expr_colon <- exprs(gse2)
pheno_colon <- pData(gse2)
expr_colon <- subset_matrix(expr_colon, n_genes = 30, n_samples = 10)
pheno_colon <- pheno_colon[colnames(expr_colon), , drop = FALSE]
save(expr_colon, pheno_colon, file = "expr_colon.rda")

# 3. GSE40595 (ovary)
gse3 <- getGEO("GSE40595", GSEMatrix = TRUE)[[1]]
expr_ovary <- exprs(gse3)
pheno_ovary <- pData(gse3)
expr_ovary <- subset_matrix(expr_ovary, n_genes = 30, n_samples = 10)
pheno_ovary <- pheno_ovary[colnames(expr_ovary), , drop = FALSE]
save(expr_ovary, pheno_ovary, file = "expr_ovary.rda")

# 4. E-MTAB-8632 (co-culture)
#ae <- ArrayExpress("E-MTAB-8632")
df <- read.delim("norm_e_mtab_8632.txt", check.names = FALSE)
df_unique <- df[!duplicated(df[,1]), ]
rownames(df_unique) <- df_unique[,1]
expr_coculture <- df_unique[,-1]
expr_coculture <- as.matrix(expr_coculture)
mode(expr_coculture) <- "numeric"
pheno_coculture <- read.delim("E-MTAB-8632.sdrf.txt", row.names = 1)
expr_coculture <- subset_matrix(expr_coculture, n_genes = 30, n_samples = 10)
pheno_coculture <- pheno_coculture[colnames(expr_coculture), , drop = FALSE]
save(expr_coculture, pheno_coculture, file = "expr_coculture.rda")

# Sposta i file .rda nella cartella data/ del pacchetto
file.copy("expr_breast.rda", "../data/expr_breast.rda", overwrite = TRUE)
file.copy("pheno_breast.rda", "../data/pheno_breast.rda", overwrite = TRUE)
file.copy("expr_colon.rda", "../data/expr_colon.rda", overwrite = TRUE)
file.copy("pheno_colon.rda", "../data/pheno_colon.rda", overwrite = TRUE)
file.copy("expr_ovary.rda", "../data/expr_ovary.rda", overwrite = TRUE)
file.copy("pheno_ovary.rda", "../data/pheno_ovary.rda", overwrite = TRUE)
file.copy("expr_coculture.rda", "../data/expr_coculture.rda", overwrite = TRUE)
file.copy("pheno_coculture.rda", "../data/pheno_coculture.rda", overwrite = TRUE)

cat("Subset di dataset tematici creati e copiati in data/.\n")
