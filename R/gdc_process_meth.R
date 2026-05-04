## max. code width ============================================================
## Libraries ==================================================================
suppressPackageStartupMessages(
  {
    library(CrossAncestryGenPhen)
  }
)

## Configs ====================================================================
configs <- list.files(
  file.path("configs", "data"), 
  pattern = "config_meth.json$",
  recursive  = TRUE, 
  full.names = TRUE
)

## Loop: All studies ==========================================================
for (config_file in configs) {

  # Read config
  cfg <- jsonlite::read_json(config_file)

  # Output directory
  if (!dir.exists(cfg$out_dir)) dir.create(cfg$out_dir, recursive = TRUE)

  # Download manifest from GDC (gdc-client)
  file_map <- gdc_download_manifest(
    manifest_file = cfg$GDC$manifest_file,
    metadata_file = cfg$GDC$metadata_file,
    client_path = "/usr/local/bin/gdc-client",
    file_dir = cfg$GDC$file_dir,
    download = cfg$GDC$download,
    verbose = FALSE
  )
  # Download clinical data from cBioportal
  clinical_data <- cbioportal_download_clinical_data(
    study_id = cfg$cBio$study_id,
    base_url = "https://www.cbioportal.org/api",
    verbose = FALSE
  )

  # Merge clinical data
  file_map <- gdc_add_clinical(
    file_map = file_map,
    clinical_data = clinical_data,
    vial = cfg$cBio$vial,
    plate = cfg$cBio$plate,
    verbose = FALSE
  )

  # Download ancestry from https://doi.org/10.1016/j.ccell.2020.04.012
  ancestry_data <- gdc_download_ancestry(
    url = "https://ars.els-cdn.com/content/image/1-s2.0-S1535610820302117-mmc2.xlsx",
    sheet = "S1 Calls per Patient",
    verbose = FALSE
  )

  # Merge ancestry data
  file_map <- gdc_add_ancestry(
    file_map = file_map,
    ancestry_data = ancestry_data,
    verbose = FALSE
  )

  # Map the meta data
  file_map <- gdc_map_metadata(
    file_map = file_map,
    map = read.csv(cfg$meta_map$map),
    verbose = FALSE
  )

  # Save metadata
  saveRDS(
    file_map, 
    file = file.path(
      cfg$out_dir, 
      "metadata_mrna.rds"
    )
  )

  # Plot
  p <- plot_gdc_demographics(
    file_map = file_map,
    attributes = cfg$demographics$attributes,
    facetting = "POOLED_GENETIC_ANCESTRY"
  )

  # Save
  ggsaveDK(
    plot = p,
    file = file.path(
      cfg$out_dir, 
      "demo_mrna.png"
    ),
    trimmed = FALSE,
    height = 10,
    width  = 20
  )

  # Beta matrix
  full_betas <- gdc_build_beta_matrix(
    file_map = file_map,
    file_dir = cfg$GDC$file_dir,
    verbose = FALSE
  )

  # Save
  saveRDS(
    full_betas, 
    file = file.path(
      cfg$out_dir, 
      "full_matrix_meth.rds"
    )
  )

  # cBioportal probes
  probe2gene <- cbioportal_download_450kprobes(
    study_id = config$cBio$study_id,
    profile_id = "data_methylation_hm27_hm450_merged.txt",
    base_url = "https://cbioportal-datahub.s3.amazonaws.com",
    verbose = FALSE
  )

  # Save
  write.csv(
    probe2gene, 
    file = file.path(
      cfg$out_dir, 
      "probe2gene.csv"
    ), 
    row.names = FALSE
  )


  # Gene-level probes
  gene_level_betas <- gdc_gene_level_beta_matrix(
    meth = full_betas,
    probe2gene = probe2gene,
    verbose = TRUE
  )

  # Save
  saveRDS(
    gene_level_betas, 
    file = file.path(
      cfg$out_dir, 
      "gene_level_matrix_meth.rds"
    )
  )
}