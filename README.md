# tkniModel

A specialized toolkit for neuroimaging analysis on local and HPC environments.

## The Workflow
1. **Voxel-wise Modeling (`modelVoxel`)**: Connects BIDS-formatted NIfTI data to an R data frame, manages fast I/O, and runs parallel models (LM, LMER, etc.) across the brain.
2. **Clustering & Thresholding (`modelCluster`)**: Identifies significant clusters based on voxel-wise and peak thresholds.
3. **Advanced Visualization (`drawMontage`)**: Creates publication-quality figures with multi-layer stacking and automated slice selection.
4. **Data Extraction (`getNIIData`)**: Extracts raw participant data from clusters or coordinates for post-hoc reporting.

## Core Functions
* **`modelVoxel()`**: Parallel statistical engine.
* **`modelCluster()`**: Comprehensive thresholding pipeline.
* **`getNIIData()`**: Fast in-memory ROI/voxel data extraction.
* **`drawMontage()`**: Grid assembly for publication figures.
* **`overlayPNG()`**: Multi-layer image compositor.
* **`slicePNG()`**: Core rendering engine.

***

## Installation
```r
devtools::install_github("TKoscik/tkniModel")
```

***

## How to Run
Since RStudio has substantial process overhead, we reccommend running modelVoxel.R using base R tools such as the linux terminal.  
Likewise the number of workers that can be effectively used depends on the number of processes R is allowed to start.  
We recommend you start an R session this way:  
```bash
R --max-connections=1024
```
Or in a script:  
```bash
Rscript --max-connections=1024 my_script.R
```
Note, that if workers are killed by your system, the overall process will check the log.nii file and re-initialize workers. This will however slow the overall processing down.  

## Example Workflow

```r
rm(list=ls())
gc()

library(tkniModel)

# Setup
dir_project <- "~/my_study"
dir_analysis <- "~/my_study/analysis/v01"

# 1. Voxel-wise Modeling
modelVoxel(
  nii_data = paste0(dir_project, "/derivatives/nifti:measure-NDI.nii.gz"),
  df_data = paste0(dir_project, "/clinical_data.csv"),
  id_var = "participant_id",
  var_factor = "group:Control,Patient",
  model_pfx = "model-NDI",
  model_fcn = function(df) {
    mdl <- lm(nii ~ age + group, data = df)
    return(as.data.frame(summary(mdl)$coefficients))
  },
  roi = paste0(dir_project, "/templates/brain_mask.nii.gz"),
  dir_save = dir_analysis,
  dir_scratch = "/tmp/tkni_scratch"
)

# 2. Clustering
modelCluster(
  nii_estimate = paste0(dir_analysis, "/model-NDI/modelResult_Estimate.nii.gz"),
  nii_test     = paste0(dir_analysis, "/model-NDI/modelResult_tvalue.nii.gz"),
  nii_pval     = paste0(dir_analysis, "/model-NDI/modelResult_Prt.nii.gz"),
  nii_mask     = paste0(dir_project, "/templates/brain_mask.nii.gz"),
  effect_name = "group", effect_volume = 3,
  peak_thresh = 0.001, cluster_size = 25,
  dir_scratch = "/tmp/tkni_cluster"
)

# 3. Visualization
drawMontage(
  bg_nii = paste0(dir_project, "/templates/T1w.nii.gz"),
  fg_nii_list = rep(paste0(dir_analysis, "/model-NDI/modelResult_tvalue.nii.gz"), 2),
  fg_mask_list = c(paste0(dir_analysis, "/model-NDI/group_neg_mask.nii.gz"),
                   paste0(dir_analysis, "/model-NDI/group_pos_mask.nii.gz")),
  fg_vol_list = c(3, 3),
  fg_colors = list(c("cyan", "blue"), c("red", "yellow")),
  layout = "9:y;9:y;9:y",
  draw_side = TRUE, draw_coords = TRUE, draw_scale = TRUE,
  file_name = "group_results_coronal",
  dir_save = dir_analysis
)

# 4. Data Extraction
cluster_data <- getNIIData(
  participant_ids = paste0("sub-", df$participant_id),
  nii_data = paste0(dir_project, "/derivatives/nifti:measure-NDI.nii.gz"),
  mask_nii = paste0(dir_analysis, "/model-NDI/group_neg_clusters.nii.gz"),
  cleanup = TRUE
)
