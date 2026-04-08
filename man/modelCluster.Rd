\name{modelCluster}
\alias{modelCluster}
\title{Neuroimaging Cluster Thresholding and Extraction}
\description{
  Performs cluster-based thresholding on neuroimaging statistical maps. The function
  identifies clusters based on a voxel-wise p-value threshold, filters them by
  spatial extent and peak significance, and outputs labeled NIfTI maps and summary tables.
}
\usage{
  modelCluster(nii_estimate, nii_test, nii_pval, nii_mask = NULL,
               effect_name = "effect", effect_volume = 1, mask_volume = 1,
               roi_thresh = 0.05, peak_thresh = 0.001, cluster_size = 25,
               connectivity = 26L, save_clusters = TRUE, save_table = TRUE,
               save_mask = TRUE, do_pos = TRUE, do_neg = TRUE,
               dir_save = getwd(), dir_scratch = NULL, cleanup = TRUE,
               verbose = TRUE)
}
\arguments{
  \item{nii_estimate}{Character. Path to the NIfTI file containing coefficient estimates (e.g., betas).}
  \item{nii_test}{Character. Path to the NIfTI file containing test statistics (e.g., T or F scores).}
  \item{nii_pval}{Character. Path to the NIfTI file containing voxel-wise p-values.}
  \item{nii_mask}{Character. Optional path to a NIfTI ROI mask.}
  \item{effect_name}{Character. Label for the effect, used in output filenames.}
  \item{effect_volume}{Integer. The volume index to read from the NIfTI files (Default: 1).}
  \item{roi_thresh}{Numeric. Initial p-value threshold for voxel inclusion (Default: 0.05).}
  \item{peak_thresh}{Numeric. P-value threshold required for at least one voxel in the cluster (Default: 0.001).}
  \item{cluster_size}{Integer. Minimum number of contiguous voxels to form a cluster (Default: 25).}
  \item{connectivity}{Integer. Neighborhood connectivity (6, 18, or 26; Default: 26).}
  \item{save_clusters}{Logical. Whether to save labeled cluster NIfTI files.}
  \item{save_table}{Logical. Whether to save a CSV summary of cluster statistics.}
  \item{save_mask}{Logical. Whether to save a binary mask of surviving clusters.}
  \item{do_pos}{Logical. Process positive effects (estimates > 0).}
  \item{do_neg}{Logical. Process negative effects (estimates < 0).}
  \item{dir_save}{Character. Directory where final outputs are saved.}
  \item{dir_scratch}{Character. Temporary directory for decompression and processing.}
  \item{cleanup}{Logical. Remove scratch files upon completion.}
  \item{verbose}{Logical. Print progress messages to the console.}
}
\details{
  The function follows a multi-step pipeline: (1) Masking by estimate sign and p-value,
  (2) Cluster identification via \code{cluster3D}, (3) Size-based filtering, (4) Peak
  coordinate and value extraction, (5) Peak-significance filtering, and (6) Renumbering
  clusters by size before final output.
}
\value{
  Invisibly returns the path to the \code{dir_save} directory.
}
\examples{
  \dontrun{
    modelCluster(
      nii_estimate = "beta_map.nii.gz",
      nii_test     = "t_stat_map.nii.gz",
      nii_pval     = "p_val_map.nii.gz",
      effect_name  = "Age_Effect",
      roi_thresh   = 0.01,
      cluster_size = 50
    )
  }
}
