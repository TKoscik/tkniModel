\name{modelVoxel}
\alias{modelVoxel}
\title{Parallel Voxel-wise Modeling for Neuroimaging Data}
\description{
  A flexible, parallelized function for running voxel-wise models on single or multiple NIfTI datasets. 
  Designed for HPC and local environments, it manages data transfer to scratch storage, 
  handles BIDS-style file matching across multiple image modalities, and executes user-defined R modeling code.
}
\usage{
modelVoxel(nii_data, df_data, id_var = "participant_id",
           var_factor = NULL, model_pfx, model_fcn, roi_nii = NULL,
           dir_save, dir_scratch = NULL, model_libraries = NULL,
           restart_log = TRUE, rand_order = TRUE,
           num_cores = parallel::detectCores() - 1,
           cleanup = TRUE, verbose = FALSE, debug = NA)
}
\arguments{
  \item{nii_data}{Character vector. Path(s) to directory/directories containing NIfTI files. Each element can optionally specify a suffix and a custom variable name separated by colons (e.g., \code{c("/path/to/niidata:_T1w:t1w", "/path/to/flairdir:_FLAIR:flair")}). Defaults variable names to \code{"nii"} for a single input or \code{"nii1"}, \code{"nii2"}, etc. for multiple inputs.}
  \item{df_data}{Character. Path to a .csv or .tsv file containing participant metadata.}
  \item{id_var}{Character vector. Variables used to match rows in \code{df_data} to filenames. Can include BIDS flag mappings (e.g., \code{c("participant_id", "session_id:ses")}).}
  \item{var_factor}{Character vector. Variables to be converted to factors. Format: \code{"var_name:levels:labels"} or \code{"var_name:levels"}.}
  \item{model_pfx}{Character. A prefix used for naming output NIfTI files and result directories.}
  \item{model_fcn}{Function. A user-defined function that takes a data frame (containing voxel values in columns named after your variables, e.g., `t1w` or `nii`) and returns a 1-row data frame of results.}
  \item{roi_nii}{Character. Optional path to a NIfTI ROI mask.}
  \item{dir_save}{Character. The permanent directory where final results will be stored.}
  \item{dir_scratch}{Character. Optional directory for fast I/O during processing. Defaults to a sub-folder in \code{dir_save}.}
  \item{model_libraries}{Character vector. Additional R packages required by \code{model_fcn}.}
  \item{restart_log}{Logical. If TRUE, restarts the nifti.io processing log.}
  \item{rand_order}{Logical. If TRUE, randomizes the order of voxel processing for better load balancing.}
  \item{num_cores}{Integer. Number of CPU cores to use for parallel processing.}
  \item{cleanup}{Logical. If TRUE, deletes the scratch directory after processing.}
  \item{verbose}{Logical. If TRUE, prints progress messages.}
  \item{debug}{Integer/NA. If an integer is provided, the function runs that many voxels sequentially for debugging.}
}
\details{
  The function copies necessary files across all specified datasets to a scratch directory to minimize network overhead. 
  It automatically handles .nii.gz decompression into the scratch space and populates columns in the participant dataframe 
  corresponding to each loaded dataset for access inside \code{model_fcn}.
}
\value{
  Returns the character path to the final results directory.
}
\examples{
\dontrun{
# Example usage with multiple datasets and custom variable names:
modelVoxel(
  nii_data = c("data/t1_dir:_T1w:t1w", "data/flair_dir:_FLAIR:flair"),
  df_data = "participants.tsv",
  id_var = c("participant_id", "session_id:ses"),
  model_pfx = "MultiModalEffect",
  dir_save = "results/",
  model_fcn = function(df) {
    fit <- lm(t1w ~ flair + age + sex, data = df)
    stats <- summary(fit)$coefficients["flair", ]
    data.frame(t_val = stats["t value"], p_val = stats["Pr(>|t|)"])
  }
)
}
}
