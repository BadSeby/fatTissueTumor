# data-raw/prepare_tcga_subset.R

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")

library(TCGAbiolinks)
library(SummarizedExperiment)

# Parametri
project <- "TCGA-BRCA" # Cambia con "TCGA-COAD" o altro se vuoi
n_genes <- 30
n_samples <- 10

# 1. Query e download (può richiedere molto tempo e banda!)
query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)
GDCdownload(query)
se <- GDCprepare(query)

# 2. Estrai matrice di espressione e fenotipo
expr_tcga <- assay(se)
pheno_tcga <- as.data.frame(colData(se))

# 3. Crea subset casuale
set.seed(123)
genes <- sample(rownames(expr_tcga), min(n_genes, nrow(expr_tcga)))
samples <- sample(colnames(expr_tcga), min(n_samples, ncol(expr_tcga)))
expr_tcga_subset <- expr_tcga[genes, samples, drop = FALSE]
pheno_tcga_subset <- pheno_tcga[samples, , drop = FALSE]

# 4. Salva i file .rda
save(expr_tcga_subset, pheno_tcga_subset, file = "expr_tcga_brca.rda")

# 5. (Opzionale) Copia in data/
file.copy("expr_tcga_brca.rda", "../data/expr_tcga_brca.rda", overwrite = TRUE)

cat("Subset TCGA creato e copiato in data/.\n")
