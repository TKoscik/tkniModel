\name{prepNII}
\alias{prepNII}
\title{Decompress and Prepare NIfTI Files for Processing}
\description{
  A utility function that handles the decompression of .nii.gz files and copies
  them to a local scratch directory. This ensures that the main rendering
  functions always work with uncompressed .nii files, which significantly
  improves I/O performance on high-performance computing (HPC) clusters.
}
\usage{
  prepNII(path, tag, dir_scratch)
}
\arguments{
  \item{path}{Character. The path to the source NIfTI file (either .nii or .nii.gz).}
  \item{tag}{Character. A short prefix (e.g., "bg", "fg1") prepended to the
    filename in the scratch directory to prevent file name collisions when
    multiple layers share the same original filename (e.g., "mask.nii.gz").}
  \item{dir_scratch}{Character. The destination directory where the prepared
    file will be stored.}
}
\details{
  The function performs the following checks:
    \enumerate{
      \item Returns \code{NULL} if the path is \code{NULL}, empty, or the string \code{"none"}.
      \item Checks if the destination file already exists in the scratch directory
      to avoid redundant decompression/copying.
      \item Uses \code{R.utils::gunzip} to decompress \code{.gz} files while
      keeping the original source file intact (\code{remove = FALSE}).
      \item Copies uncompressed \code{.nii} files directly to the scratch
      location if they are not already there.
    }
}
\value{
  Character. The full path to the prepared .nii file in the scratch directory,
  or \code{NULL} if no valid input path was provided.
}
\note{
  This is a helper function primarily used internally by \code{\link{slicePNG}}
  and \code{\link{overlayPNG}}.
}
