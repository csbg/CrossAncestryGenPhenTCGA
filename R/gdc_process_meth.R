library(CrossAncestryGenPhen)

# ===== Load config =====
config <- jsonlite::read_json(
  "download/GDCportal/TCGA_THCA/config_meth.json", 
  simplifyVector = TRUE
)


# ===== Download from GDCportal =====
## --- Download manifest ---
file_map <- gdc_download_manifest(
  manifest_file = config$GDC$manifest_file,
  metadata_file = config$GDC$metadata_file,
  client_path   = "/usr/local/bin/gdc-client",
  file_dir      = config$GDC$file_dir,
  download      = config$GDC$download,
  verbose       = TRUE
)


# ===== Create meta data =====
## --- Donwload clinical data ---
clinical_data <- cbioportal_download_clinical_data(
  study_id = config$cBio$study_id,
  base_url = "https://www.cbioportal.org/api",
  verbose  = TRUE
)


## --- Merge with clinical data ---
file_map <- gdc_add_clinical(
  file_map      = file_map,
  clinical_data = clinical_data,
  vial          = config$cBio$vial,
  plate         = config$cBio$plate,
  verbose       = TRUE
)


## --- Download ancestry call ---
ancestry_data <- gdc_download_ancestry(
  url     = "https://ars.els-cdn.com/content/image/1-s2.0-S1535610820302117-mmc2.xlsx",
  sheet   = "S1 Calls per Patient",
  verbose = TRUE
)


## --- Add ancestry ---
file_map <- gdc_add_ancestry(
  file_map      = file_map,
  ancestry_data = ancestry_data,
  verbose       = TRUE
)


## --- Map metadata values ---
file_map <- gdc_map_metadata(
  file_map = file_map,
  map      = read.csv(config$meta_map$map),
  verbose  = TRUE
)
saveRDS(file_map, file = file.path(config$out_dir, "metadata_meth.rds"))


## --- Demographics ---
p <- plot_gdc_demographics(
  file_map   = file_map,
  attributes = config$demographics$attributes,
  facetting  = "POOLED_GENETIC_ANCESTRY"
)
ggsaveDK(
  plot = p,
  file = file.path(config$out_dir,"demo_meth.png"),
  trimmed = FALSE,
  height = 10,
  width  = 20
)


## --- Create beta matrix ---
full_betas <- gdc_build_beta_matrix(
  file_map = file_map,
  file_dir = config$GDC$file_dir,
  verbose  = TRUE
)
saveRDS(full_betas, file = file.path(config$out_dir, "full_matrix_meth.rds"))


## --- cBioportal probes ---
probe2gene <- cbioportal_download_450kprobes(
  study_id = config$cBio$study_id,
  profile_id = "data_methylation_hm27_hm450_merged.txt",
  base_url = "https://cbioportal-datahub.s3.amazonaws.com",
  verbose = TRUE
)
write.csv(probe2gene, file = file.path(config$out_dir, "probe2gene.csv"), row.names = FALSE)


## --- Gene-level probes ---
gene_level_betas <- gdc_gene_level_beta_matrix(
  meth = full_betas,
  probe2gene = probe2gene,
  verbose = TRUE
)
saveRDS(gene_level_betas, file = file.path(config$out_dir, "gene_level_matrix_meth.rds"))
