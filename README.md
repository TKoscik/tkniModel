# tkniModel

`tkniModel` is an R package designed for high-performance neuroimaging analysis. It is specifically optimized for Cluster/HPC environments where local scratch storage and parallel processing are essential for scaling voxel-wise models.

## The Workflow

The package is designed to be used in three distinct phases:

1.  **Voxel-wise Modeling (`modelVoxel`)**: Connects your BIDS-formatted NIfTI data to an R data frame, manages fast I/O via scratch directories, and runs parallel models (LM, LMER, etc.) across the brain.
2.  **Clustering & Extraction (`modelCluster`)**: Identifies significant clusters based on voxel-wise and peak thresholds, sorts them by size, and generates summary statistics and labeled NIfTI maps.
3.  **Visualization (`drawOverlay` & `drawEffect`)**: Produces publication-quality figures, including multi-layered brain montages with real-world coordinates and adjusted effect plots that account for random-effects variance.

---

## Core Functions

### Analysis & Processing
*   **`modelVoxel()`**: The primary engine. It maps participant metadata to NIfTI files, handles local decompression on compute nodes, and executes user-defined R functions across every voxel.
*   **`modelCluster()`**: A comprehensive thresholding pipeline. It identifies contiguous clusters, applies size and peak-significance filters, and renumbers clusters for easy reporting.
*   **`cluster3D()`**: A memory-efficient 3D seed-fill algorithm for identifying connected components using 6, 18, or 26-neighbor connectivity.
*   **`clusterTable()`**: Extracts peak intensities, world coordinates, and center-of-gravity stats for labeled clusters.

### Visualization
*   **`drawOverlay()`**: Generates "lightbox" montage figures. Automatically finds slices with data, stacks multiple statistical layers over anatomy, and calculates coordinate labels from NIfTI `sform` headers.
*   **`drawEffect()`**: Plots the relationship between behavioral variables and brain data. For mixed-effects models, it "partials out" random effects to show the true fixed-effect relationship.
*   **`slicesToRaster()`**: A utility that transforms 3D volumes into 2D montage grids compatible with `ggplot2`.

---

## Installation

```r
# Install from GitHub
devtools::install_github("TKoscik/tkniModel")
```

## Workflow Example
```r
library(tkniModel)

# --- 1. Project Setup ---
# Define paths to your data and analysis folders
dir_project  <- "~/my_study"
dir_analysis <- "~/my_study/analysis/v01"
dir_scratch  <- "/tmp/tkni_scratch"

# --- 2. Run Voxel-wise Modeling ---
# This example tests for 'age' and 'group' effects across all brain voxels
modelVoxel(
  # Path to participant NIfTIs (using path:suffix notation)
  nii_data = paste0(dir_project, "/derivatives/nifti:measure-NDI.nii.gz"),
  
  # Path to your clinical/behavioral CSV
  df_data = paste0(dir_project, "/clinical_data.csv"),
  
  # Specify the participant ID column and any factor levels
  id_var = "participant_id",
  var_factor = "group:Control,Patient",
  
  # Define the statistical model to run at each voxel
  model_pfx = "model-NDI_age_by_group",
  model_fcn = function(df) {
    # 'nii' is the reserved variable name for the voxel intensity
    mdl <- lm(nii ~ age + group, data = df)
    
    # Return the coefficients as a data frame
    stats <- as.data.frame(summary(mdl)$coefficients)
    return(stats)
  },
  
  # Optional: Restrict analysis to a specific ROI or brain mask
  roi = paste0(dir_project, "/templates/brain_mask.nii.gz"),
  
  # Directory management
  dir_save = dir_analysis,
  dir_scratch = dir_scratch
)

# --- 3. Thresholding and Clustering ---
# Convert raw voxel results into significant clusters 
# based on p-value and cluster size (extent) thresholds.
modelCluster(
  # Paths to the NIfTI outputs from modelVoxel
  nii_estimate = paste0(dir_analysis, "/model-NDI/modelResult_Estimate.nii.gz"),
  nii_test     = paste0(dir_analysis, "/model-NDI/modelResult_tvalue.nii.gz"),
  nii_pval     = paste0(dir_analysis, "/model-NDI/modelResult_Prt.nii.gz"),
  nii_mask     = paste0(dir_project, "/templates/brain_mask.nii.gz"),
  
  # Select the specific effect from your model (e.g., the 3rd volume/coefficient)
  effect_name   = "group",
  effect_volume = 3,
  
  # Set thresholds: 
  # Clusters must have a peak p < 0.001 and at least 25 voxels
  roi_thresh   = 0.05,
  peak_thresh  = 0.001,
  cluster_size = 25,
  connectivity = 26, # 26-neighbor connectivity (faces, edges, & corners)
  
  # Search for both positive and negative effects
  do_pos = TRUE, 
  do_neg = TRUE,
  
  dir_scratch = "/tmp/tkni_cluster_temp"
)

# --- 4. Generate Publication-Ready Montage ---
# Create a multi-slice grid showing significant clusters 
# overlaid on anatomical background with ROI outlines.
drawMontage(
  # Background Anatomy
  bg_nii = paste0(dir_project, "/templates/T1w_template.nii.gz"),
  bg_mask = paste0(dir_project, "/templates/brain_mask.nii.gz"),
  bg_threshold_pct = c(0.025, 0.975),
  
  # Foreground Layers (Statistical Maps)
  # Here we overlay the same t-value map twice, using different masks 
  # to separate and color-code negative vs. positive effects.
  fg_nii_list = rep(paste0(dir_analysis, "/model-NDI/modelResult_tvalue.nii.gz"), 2),
  fg_mask_list = c(
    paste0(dir_analysis, "/model-NDI/effect-group_dir-neg_mask.nii.gz"),
    paste0(dir_analysis, "/model-NDI/effect-group_dir-pos_mask.nii.gz")
  ),
  fg_vol_list = c(3, 3),        # Volume index for the 'group' effect
  fg_alphas   = 0.5,            # Set overlay transparency
  fg_colors   = list(c("cyan", "blue"), c("red", "yellow")),
  
  # ROI Outlines (Pink for negative clusters, Green for positive)
  roi_nii_list = c(
    paste0(dir_analysis, "/model-NDI/effect-group_dir-neg_mask.nii.gz"),
    paste0(dir_analysis, "/model-NDI/effect-group_dir-pos_mask.nii.gz")
  ),
  roi_color = c("hotpink", "limegreen"),
  roi_value = "all",
  
  # Layout and Annotations
  layout    = "9:y;9:y;9:y",    # 3 rows of 9 coronal slices each
  edge_clip = 0,               # Don't clip edges of the bounding box
  draw_side   = TRUE,          # Add laterality arrows (<- L)
  draw_coords = TRUE,          # Add world coordinates (mm) to each slice
  draw_scale  = TRUE,          # Add physical scale bar
  
  # Color Bar configuration
  draw_cbar     = "horizontal",
  cbar_location = "southeast",
  
  # File Management
  file_name   = "effect-group_significant_results",
  dir_save    = dir_analysis,
  dir_scratch = "/tmp/montage_temp",
  cleanup     = TRUE           # Automatically wipe temp files after finishing
)

# --- 5. Extract Cluster Data for Post-hoc Analysis ---
# Pull the mean values from each significant cluster for 
# every participant to create scatter plots or boxplots.

# Load your original behavioral/clinical dataframe
df <- read.csv(paste0(dir_project, "/clinical_data.csv"))

# Extract brain data from the cluster mask generated in Step 3
# getNIIData will find the correct files and calculate means in RAM
cluster_data <- getNIIData(
    # Provide the BIDS identifiers
    participant_ids = paste0("sub-", df$participant_id),
    
    # Path to participant NIfTIs (using path:suffix notation)
    nii_data = paste0(dir_project, "/derivatives/nifti:measure-NDI.nii.gz"),
    nii_vol  = 1,
    
    # The mask containing labeled clusters (e.g., cluster 1, 2, 3...)
    mask_nii   = paste0(dir_analysis, "/model-NDI/effect-group_neg_clusters.nii.gz"),
    mask_vol   = 1,
    mask_value = NULL, # Extract all labeled clusters found in the mask
    
    cleanup = TRUE,
    verbose = TRUE
)

# Merge the brain data back into your clinical dataframe for final plotting
final_df <- merge(df, cluster_data, by.x = "participant_id", by.y = "participant_id")

# Now you're ready for ggplot2!
# ggplot(final_df, aes(x=group, y=cluster_001)) + geom_boxplot()


```
