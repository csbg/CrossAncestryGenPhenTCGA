# Portability of cancer effects across genetic ancestries
This github hosts the code for the publication: **insert link**

## Reporting summary
Additional information for the reporting can be found in the directory **paper**.
```
paper/ 
├── reporting_summary.txt
│
├── summary_subset_interaction_effect/ 
│   └── summary_sample_n.csv             # Number of samples per comaprison
│
├── summary_subset_prediction_effect/
│   └── summary_sample_n.csv  
│
└── summary_limma_interaction_effect/
    └── summary_sample_n.csv
```
Same samples are used across all analysis. For subsetting analysis (subset_interaction_effect_script.R & subset_prediction_effect_script.R) the same subsets (same patient ids) are reused!


## Data download
The molecular data to run this script is downloaded from [GDC Data Portal](https://portal.gdc.cancer.gov) using the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) with the **manifest.txt**, **metadata.json** and the corresponding **config_<molecular_level>.json** files in *configs/data/*. 
```
configs/ 
└── data/ 
    ├── TCGA_BRCA/ 
    │    ├── GDC_manifest_TCGA_BRCA_meth.txt     
    │    ├── GDC_manifest_TCGA_BRCA_mrna.txt     
    │    ├── GDC_metadata_TCGA_BRCA_meth.json    
    │    ├── GDC_metadata_TCGA_BRCA_mrna.json    
    │    ├── config_meth.json            # Config for script `gdc_process_meth.R`
    │    ├── config_mrna.json            # Config for script `gdc_process_mrna.R`
    │    └── metadata_map.csv            # Map to rename column names and values
    │
    ├── TCGA_THCA/ 
    │    └── ...     
    │
    └── TCGA_UCEC/ 
         └── ...
```
Once you have installed the `gdc-client` you can run the scripts `gdc_process_meth.R`, `gdc_process_meth.R`. Please be aware that the scripts ask for the location of the client. The scripts will download the molecular data (from GDC Portal), additional clinical data (cBioportal) and genetic ancestry (https://doi.org/10.1016/j.ccell.2020.04.012). Once the data for each molecular layer (methylation & expression) is downloaded the data is combined into one object using `gdc_create_study.R`. The scripts will create a new directory *download/*.
```
download/ 
└── GDCportal/ 
    ├── TCGA_BRCA/ 
    │    ├── meth_files/                 # Raw patient files
    │    ├── mrna_files/     
    │    ├── demo_meth.png               # Demographic plots
    │    ├── demo_mrna.png       
    │    ├── full_matrix_meth.rds        # All downloaded features                    
    │    ├── full_matrix_mrna.rds
    │    ├── gene_level_matrix_meth.rds  # Summarized on gene level
    │    ├── gene_level_matrix_mrna.rds    
    │    ├── metadata_meth.rds           # Metadata files
    │    ├── metadata_mrna.rds
    │    ├── probe2gene.csv              # Mapping of probe to gene
    │    └── TCGA_BRCA_omics_layers.rds  # Combined study object (meth & mrna)
    │
    ├── TCGA_THCA/ 
    │    └── ...
    │
    └── TCGA_UCEC/ 
         └── ...
```
## Data simulation
Both data simulations (synthetic and permuted) were performed using data from TCGA BRCA. The configs for these analysis are in *configs/sims/*.
The scripts `sim_synthetic_script.R` and `sim_permutation_script.R` will use the same config files and create *results/sims/analysis*. The scripts `sim_synthetic_figures.R` and `sim_permutation_figures.R` will create the figures for the publication.
```
results/ 
└── sims/ 
    ├── analysis/ 
    │    ├── TCGA_BRCA/
    │    │    ├── permutation_sim/
    │    │    │    ├── res_file1            # One res file per config
    │    │    │    └── ...
    │    │    │
    │    │    └── synthetic_sim/
    │    │         ├── res_file1            # One res file per config
    │    │         └── ...
    │    │
    │    └── TCGA_UCEC/ 
    │         └── ...
    │
    ├── figures/ 
    └── tables/
```
## Differential expression and ML prediction
Both analysis differntial expression and ML prediction use the same config files in *conigs/tcga/*.
The scripts `subset_interaction_effect_script.R` and `subset_prediction_effect_script.R` will use the same config files and create *results/tcga/analysis*. 
The scripts `subset_interaction_effect_figures.R` and `subset_prediction_effect_figures.R`will create the figures for the publication.
```
results/ 
└── tcga/ 
    ├── analysis/ 
    │    ├── TCGA_BRCA/
    │    │    ├── subset_interaction_effect/
    │    │    │    ├── res_file1            # One res file per config
    │    │    │    └── ...
    │    │    │
    │    │    └── subset_prediction_effect/
    │    │         ├── res_file1            # One res file per config
    │    │         └── ...
    │    │
    │    ├── TCGA_THCA/ 
    │    │    └── ...
    │    │
    │    └── TCGA_UCEC/ 
    │          └── ...
    │
    ├── figures/ 
    └── tables/
```
## Requirements
Before running these scripts, please install the **CrossAncestryGenPhen** package from GitHub.

Link: https://github.com/DKatzlberger/CrossAncestryGenPhen

You can do this in R using:

```r
install.packages("devtools")  # if not already installed
devtools::install_github("DKatzlberger/CrossAncestryGenPhen")
```

Make sure the installation completes successfully, then load the package with:

```r
library(CrossAncestryGenPhen)
```

Once the package is installed and loaded, the scripts should run as intended.
