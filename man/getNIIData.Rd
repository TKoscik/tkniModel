\name{getNIIData}
\alias{getNIIData}
\title{Optimized Extraction of Brain Data by Participant ID}
\description{
  Finds participant NIfTI files based on BIDS identifier strings and extracts
  brain data. The function supports extracting mean values from labeled clusters
  within a mask or specific voxel values from a coordinate matrix. It is
  optimized for speed by loading the entire participant volume into memory
  once per subject and applying ROI masks or coordinate lookups in RAM.
}
\usage{
getNIIData(participant_ids, nii_data, mask_nii = NULL, mask_vol = 1,
           mask_value = NULL, coords = NULL, nii_vol = 1,
           dir_scratch = NULL, cleanup = TRUE, verbose = TRUE)
}
\arguments{
  \item{participant_ids}{Character vector. Full BIDS identifier strings for
    participants (e.g., \code{c("sub-101_ses-1", "sub-102_ses-1")}).}
  \item{nii_data}{Character. The path to the data directory or a search string
    combining path and suffix (e.g., \code{"path/to/data:_T1w.nii.gz"}).}
  \item{mask_nii}{Character. Optional path to a labeled cluster mask file.}
  \item{mask_vol}{Integer. The volume index to use if the mask file is
    multi-volume (defaults to 1).}
  \item{mask_value}{Numeric vector. Specific label values in the mask to
    extract. If \code{NULL}, all non-zero labels are extracted.}
  \item{coords}{Matrix. An \eqn{N \times 3} matrix of voxel coordinates for
    specific voxel extraction.}
  \item{nii_vol}{Integer. The volume index to extract from for 4D participant
    data (defaults to 1).}
  \item{dir_scratch}{Character. Path to a scratch directory for decompression.
    If \code{NULL}, a temporary directory is used.}
  \item{cleanup}{Logical. If \code{TRUE} (default), deletes temporary unzipped
    NIfTI files from the scratch directory after extraction is complete.
    Set to \code{FALSE} if you intend to reuse unzipped files for subsequent
    analyses.}
  \item{verbose}{Logical. Whether to print progress messages to the console.}
}
\details{
  The function uses a hierarchical approach to identify NIfTI files based on
  BIDS strings and suffixes. For performance, it utilizes \code{\link[nifti.io]{read.nii.volume}}
  to load participant data into an R array once per subject.

  When a mask is provided, each unique label (\code{mask_value}) is extracted
  as a separate column in the returned data frame, named \code{cluster_#}.
  When \code{coords} are provided, the function treats the coordinates as a
  single ROI and returns the mean value in a column named \code{cluster_coords}.
}
\value{
  A \code{data.frame} containing the original \code{participant_id}, the matched
  \code{nii_file} path, and additional columns for each extracted ROI or
  coordinate value.
}
\examples{
\dontrun{
# Extract mean NDI values from clusters 1 and 2 and cleanup scratch
results <- getNIIData(
  participant_ids = c("sub-01", "sub-02"),
  nii_data = "/path/to/data:_NDI.nii.gz",
  mask_nii = "path/to/mask.nii.gz",
  mask_value = c(1, 2),
  cleanup = TRUE
)

# Extract voxel data without cleaning up (faster for repeated calls)
vox_data <- getNIIData(
  participant_ids = c("sub-01", "sub-02"),
  nii_data = "/path/to/data:_T1w.nii.gz",
  coords = matrix(c(45, 60, 32), ncol = 3),
  cleanup = FALSE
)
}
}
