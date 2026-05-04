## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
  }
)

## Study dirs =================================================================
studies <- list.dirs(
  file.path("download", "GDCportal"),
  recursive = FALSE,
  full.names = TRUE
)

## Loop: All studies ==========================================================
for (study_dir in studies) {

  # Meth
  meth_layer <- gdc_omics_layer(
    meta = readRDS(file.path(study_dir, "metadata_meth.rds")),
    matr = readRDS(file.path(study_dir, "gene_level_matrix_meth.rds")),
    id = "SAMPLE_ID",
    tech = "meth",
    verbose = FALSE
  )

  # Expr
  mrna_layer <- gdc_omics_layer(
    meta = readRDS(file.path(study_dir, "metadata_mrna.rds")),
    matr = readRDS(file.path(study_dir, "gene_level_matrix_mrna.rds")),
    id = "SAMPLE_ID",
    tech = "mrna",
    verbose = FALSE
  )

  # Combine into study
  study_obj  <- c(meth_layer, mrna_layer)
  study_name <- basename(study_dir)

  # Save
  saveRDS(
    study_obj,
    file = file.path(
      study_dir, 
      paste0(
        study_name, 
        "_omics_layers.rds"
      )
    )
  )
}