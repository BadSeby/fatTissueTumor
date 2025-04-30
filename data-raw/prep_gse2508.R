# Scarica la serie GEO
gse <- GEOquery::getGEO("GSE2508", GSEMatrix = TRUE)[[1]]

# Estrai la matrice di espressione
expr <- exprs(gse)
print(dim(expr)) # Controlla dimensioni

# Riduci la dimensione per esempio (es. primi 100 geni e 10 campioni)
n_geni <- min(100, nrow(expr))
n_campioni <- min(10, ncol(expr))
expr_small <- expr[1:n_geni, 1:n_campioni]

# Salva il dataset in formato RData per il pacchetto
usethis::use_data(expr_small, overwrite = TRUE)
