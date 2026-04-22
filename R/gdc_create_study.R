library(CrossAncestryGenPhen)

# ===== Data =====
# --- Methylation ---
meth_layer <- gdc_omics_layer(
  meta    = readRDS("download/GDCportal/TCGA_THCA/metadata_meth.rds"),
  matr    = readRDS("download/GDCportal/TCGA_THCA/gene_level_matrix_meth.rds"),
  id      = "SAMPLE_ID",
  tech    = "meth",
  verbose = TRUE
)

# --- Expression ---
mrna_layer <- gdc_omics_layer(
  meta    = readRDS("download/GDCportal/TCGA_THCA/metadata_mrna.rds"),
  matr    = readRDS("download/GDCportal/TCGA_THCA/gene_level_matrix_mrna.rds"),
  id      = "SAMPLE_ID",
  tech    = "mrna",
  verbose = TRUE
)

# ===== Study =====
study <- c(meth_layer, mrna_layer)

# Save
saveRDS(study, file.path("download/GDCportal/TCGA_THCA", "TCGA_THCA_omics_layers.rds"))
