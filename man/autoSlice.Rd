\name{autoSlice}
\alias{autoSlice}
\title{Automated Slice Selection for Neuroimaging Montages}
\description{
  Calculates optimal slice indices for anatomical montages based on the spatial
  extent of statistical results. It parses a custom layout string and returns
  centered slice indices for each requested plane.
}
\usage{
  autoSlice(nii_paths, layout = "1:x,1:y,1:z", edge_clip = 5,
            dir_scratch = tempdir())
}
\arguments{
  \item{nii_paths}{Character vector. Paths to the NIfTI files (e.g., masks or
                                                               statistical maps) used to determine the bounding box of results.}
  \item{layout}{Character. A delimited string defining the montage structure.
    Rows are separated by \code{;}, columns by \code{,}, and slice counts per
    plane by \code{:}. Example: \code{"5:z;5:z"} creates two rows of 5 axial slices.}
  \item{edge_clip}{Integer. Number of voxels to trim from the edges of the
    detected bounding box to avoid plotting empty space or image artifacts.}
  \item{dir_scratch}{Character. Temporary directory for decompressing input files.}
}
\details{
  The function finds the union of non-zero voxels across all provided
  \code{nii_paths} to define a result-specific bounding box. It then partitions
  this box into the requested number of slices for each plane (x, y, z),
  ensuring the slices are evenly spaced and centered within the data extent.
  This mimics the automated slice selection logic found in high-performance
  bash-based neuroimaging toolkits.
}
\value{
  A named list containing numeric vectors of indices for each plane:
    \itemize{
      \item \code{x}: Sagittal slice indices
      \item \code{y}: Coronal slice indices
      \item \code{z}: Axial slice indices
    }
}
\examples{
  \dontrun{
    # Determine slices for a 3x5 axial montage based on a T-stat map
    slices <- autoSlice(
      nii_paths = "tstat_map.nii.gz",
      layout = "5:z;5:z;5:z",
      edge_clip = 10
    )
  }
}
