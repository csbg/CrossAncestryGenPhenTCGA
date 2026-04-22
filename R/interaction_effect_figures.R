library(CrossAncestryGenPhen)
library(clusterProfiler)
library(ComplexHeatmap)
library(data.table)
library(yardstick)
library(patchwork)
library(ggplot2)
library(scales)
library(edgeR)
library(Rtsne)
library(grid)

out_dir <- "figures/interaction_effect"

# ======= Color configurations =====
study_colors <- c(
  "BRCA" = "#5e3c99",   # deep indigo
  "THCA" = "#fdb863",   # peach-sand
  "UCEC" = "#b2abd2"    # light lavender
)

tech_colors <- c(
  "mrna" = "#6B8E23",  
  "meth" = "#EE9572" 
)

phenotype_colors <- c(
  "Normal vs Primary"       = "#ff7f00",   # orange
  "Basal vs non-Basal"      = "#e41a1c",   # red
  "LumA vs LumB"            = "#4daf4a",   # green
  "Serous vs Endometrioid"  = "#377eb8",   # blue
  "Classical vs Follicular" = "#984ea3",   # purple (NEW, unique)
  "M0 vs MX"                = "#a65628"    # brown-orange (NEW, unique)
)

ancestry_colors <- c(
  "EUR"   = "#0072B2",
  "AFR"   = "#D55E00",
  "EAS"   = "#56B4E9",
  "AMR"   = "#E69F00",
  "SAS"   = "#009E73",
  "ADMIX" = "#999999"
)

sample_type_colors <- c(
  "Normal"     = "#7570b3",
  "Primary"    = "#1b9e77", 
  "Metastatic" = "brown"
)

subtype_colors <- c(
  "LumA"   = "#8DB6CD",
  "LumB"   = "#e31a1c",
  "Her2"   = "#33a02c",
  "Basal"  = "#ff7f00",
  "Normal" = "#CD00CD"
)

sex_colors <- c(
  "Male"   = "#8B5742",
  "Female" = "#EEE9BF"
)

effect_group_colors <- c(
  "De novo effect"                  = "#D62728",  
  "Minor interaction effect"        = "black", 
  "non-EUR enhances\ncancer effect" = "#1F77B4",  
  "non-EUR reverts\ncancer effect"  = "#AEC7E8",  
  "Cancer enhances\nancestry effect"= "#98DF8A",  
  "Cancer reverts\nancestry effect" = "#2CA02C"
)

### ======= DGE configurations =======
results <- list(
  BRCA = "results/TCGA_BRCA",
  UCEC = "results/TCGA_UCEC",
  THCA = "results/TCGA_THCA"
)


### ====== Interaction effect ======
## --- Limma interaction results ---
interaction <- rbindlist(lapply(names(results), function(study) {
  file_path <- file.path(results[[study]], "interaction_effect")
  dge_files <- list.files(file_path, pattern = "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
  ## Load all files 
  all_files <- lapply(dge_files, readRDS)
  dt <- rbindlist(all_files, use.names = TRUE, fill = TRUE)
  ## Set phenotype
  dt[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  p <- unique(dt$phenotype)
  dt[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
  ## Return
  dt
}), fill = TRUE)

## --- Subset interaction results ---
subset_interaction <- rbindlist(lapply(names(results), function(study) {
  file_path <- file.path(results[[study]], "subset_interaction_effect")
  dge_files <- list.files(file_path, pattern = "dge_res\\.rds$", recursive = TRUE, full.names = TRUE)
  ## Load all files
  all_files <- lapply(dge_files, readRDS)
  dt <- rbindlist(lapply(all_files, \(x) x$summary_stats), use.names = TRUE, fill = TRUE)
  ## Set phenotype
  dt[, phenotype := gsub("_", "-", paste(g_1, "vs", g_2))]
  p <- unique(dt$phenotype)
  dt[, phenotype := factor(phenotype, levels = c(p[grepl("^Normal", p)], p[!grepl("^Normal", p)]))]
  ## Return
  dt
}), fill = TRUE)



### ==== Correlation with limma ====
## ---- Complex heatmap (limma-voom correaltion) ---
## Reduce dt
corr_dt <- rbind(interaction, subset_interaction, fill = TRUE)
corr_dt <- corr_dt[,p_value := -log10(p_value + 1)]
corr_dt <- corr_dt[, .(study, phenotype, method, toi, a_2, feature, p_value)]
corr_dt[, pipeline := paste(study, toi, phenotype, a_2, sep = "|")]

## Build matrix
wide <- dcast(corr_dt, pipeline + feature ~ method, value.var = "p_value")
wide_list <- split(wide, wide$pipeline)
corr_dt   <- rbindlist(lapply(names(wide_list), function(p) {
  df <- wide_list[[p]]
  m  <- as.matrix(df[, meth_cols, with = FALSE])
  cmat <- cor(m, use = "pairwise.complete.obs", method = "pearson")
  out <- as.data.table(as.table(cmat))
  setnames(out, c("Var1", "Var2", "Freq"))
  out[, pipeline := p]
  out
}))


# Filter and parse
corr_ann <- corr_dt[Var1 != Var2 & Var1 < Var2]
corr_ann[, c("study", "tech", "phenotype", "a_2") := tstrsplit(pipeline, "\\|")]

# Build metadata table indexed by pipeline
ann <- unique(corr_ann[, .(pipeline, study, toi, phenotype, a_2)])
setkey(ann, pipeline)

# Extract correlations and build matrix
heat_df <- corr_ann[, .(pipeline, corr = Freq)]
mat <- t(as.matrix(heat_df$corr))

# Set colnames to PIPELINE FIRST (this is critical!)
colnames(mat) <- heat_df$pipeline
rownames(mat) <- "subset-cct vs. limma-voom"

ann <- ann[colnames(mat)]
display_names <- ann$a_2           
names(display_names) <- ann$pipeline
colnames(mat) <- display_names


# Prepare sig_genes
sig_genes_interaction[, pipeline := paste(study, toi, phenotype, a_2, sep = "|")]
sig_vec <- sig_genes_interaction[match(heat_df$pipeline, pipeline), sig_genes]

# Barplot annotation (bars below the heatmap)
ha_bar <- HeatmapAnnotation(
  "Log10 sig. genes" = anno_barplot(
    log10(sig_vec + 1),
    gp = gpar(
      fill = "grey80",
      col  = "black"    
    ),
    bar_width = 0.6,
    border = FALSE,
    height = unit(1, "cm"),
    axis_param = list(
      gp = text_gp,
      side = "left"
    )
  ),
  annotation_name_side = "left",
  annotation_name_gp = text_gp
)

# Metadata annotation
ha_meta <- HeatmapAnnotation(
  #Ancestry  = ann$a_2,
  Phenotype = ann$phenotype,
  Tech      = ann$tech,
  Study     = ann$study,
  col = list(
    #Ancestry  = ancestry_colors[unique(ann$a_2)],
    Phenotype = phenotype_colors[unique(ann$phenotype)],
    Tech      = tech_colors[unique(ann$tech)],
    Study     = study_colors[unique(ann$study)]
  ),
  simple_anno_size = unit(8, "pt"),
  annotation_name_side = "left",
  annotation_name_gp   = text_gp,
  annotation_legend_param = list(
    title_gp  = text_gp,
    labels_gp = text_gp,
    legend_direction = "horizontal",
    nrow = 2
  )
)


ht <- Heatmap(
  mat,
  name = "Pearson",
  column_title = "Ancestry",
  column_title_side = "bottom",
  column_title_gp = text_gp,
  col  =  colorRamp2(c(0, 1), c("white", "red")),                     
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = FALSE,
  top_annotation = ha_bar,     
  bottom_annotation = ha_meta,
  height = unit(1.5, "cm"),
  width = unit(16, "cm"),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  column_names_centered = FALSE,
  heatmap_legend_param = list(
    title_gp   = text_gp,
    labels_gp  = text_gp,               
    legend_height = unit(10, "pt"),                       
    legend_width  = unit(8,  "pt"),                       
    grid_width    = unit(8,  "pt"),                       
    grid_height   = unit(8,  "pt")                        
  )
)


svg(file.path(out_dir, "heatmap.svg"), width  = 15, height = 15)
draw(
  ht, 
  show_annotation_legend = FALSE, 
  show_heatmap_legend = TRUE,
  padding = unit(c(0, 0, 0, 0), "mm")
)
dev.off()
