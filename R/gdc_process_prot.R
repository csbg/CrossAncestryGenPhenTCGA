library(CrossAncestryGenPhen)

# ===== Load config =====
config <- jsonlite::read_json("download/GDCportal/TCGA_BRCA/config_prot.json", simplifyVector = TRUE)


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


## --- Add clinical data ---
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
saveRDS(file_map, file = file.path(config$out_dir, "metadata_prot.rds"))


## --- Demographics ---
p <- plot_gdc_demographics(
  file_map   = file_map,
  attributes = config$demographics$attributes,
  facetting  = "POOLED_GENETIC_ANCESTRY"
)
ggsaveDK(
  plot = p,
  file = file.path(config$out_dir, "demo_prot.png"),
  trimmed = FALSE,
  height = 10,
  width  = 20
)


# ===== Create molecular data =====
full_rppa <- gdc_build_rppa_matrix(
  file_map = file_map,
  file_dir = "download/GDCportal/TCGA_BRCA/prot_files",
  verbose = TRUE
)
saveRDS(full_rppa, file = file.path(config$out_dir, "full_matrix_prot.rds"))
