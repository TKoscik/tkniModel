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

## Quick Example
```r
library(tkniModel)

# 1. Run the voxel-wise model
modelVoxel(
  nii_data = "path/to/niftis:_T1w",
  df_data = "participants.csv",
  id_var = c("participant_id", "session_id:ses"),
  model_pfx = "MainEffect",
  model_fcn = function(df) {
    fit <- lm(nii ~ age + sex, data = df)
    res <- summary(fit)$coefficients["age", ]
    data.frame(t_val = res["t value"], p_val = res["Pr(>|t|)"])
  },
  dir_save = "./results"
)

# 2. Cluster the results
modelCluster(
  nii_estimate = "beta.nii",
  nii_test = "tstat.nii",
  nii_pval = "p_val.nii",
  roi_thresh = 0.01,
  cluster_size = 50,
  dir_save = "./results"
)
```
