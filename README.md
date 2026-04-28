# Portability of cancer effects across genetic ancestries
This github hosts the code for the publication: **insert link**

## Data download
The molecular data to run this script is downloaded from [GDC Data Portal](https://portal.gdc.cancer.gov) using the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) with the **manifest.txt** and **metadata.json** in configs/data/:

```
configs/ 
└── data/ 
    ├── TCGA_BRCA/ 
    │       ├── GDC_manifest_TCGA_BRCA_meth.txt     # Manifest file for methylation
    │       ├── GDC_manifest_TCGA_BRCA_mrna.txt     # Manifest file for expression
    │       ├── GDC_metadata_TCGA_BRCA_meth.json    # Metadata for methylation
    │       ├── GDC_metadata_TCGA_BRCA_mrna.json    # Metadata for expression
    │       ├── config_meth.json                    # Config for script `gdc_process_mrna.R`
    │       └── config_mrna.json                    # Config for script `gdc_process_meth.R`
    ├── TCGA_THCA/ 
    │       ├── GDC_manifest_TCGA_THCA_meth.txt     
    │       ├── GDC_manifest_TCGA_THCA_mrna.txt     
    │       ├── GDC_metadata_TCGA_THCA_meth.json    
    │       ├── GDC_metadata_TCGA_THCA_mrna.json  
    │       ├── config_meth.json 
    │       └── config_mrna.json 
    └── TCGA_UCEC/ 
            ├── GDC_manifest_TCGA_UCEC_meth.txt     
            ├── GDC_manifest_TCGA_UCEC_mrna.txt     
            ├── GDC_metadata_TCGA_UCEC_meth.json    
            ├── GDC_metadata_TCGA_UCEC_mrna.json 
            ├── config_meth.json 
            └── config_mrna.json
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
