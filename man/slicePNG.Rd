\name{slicePNG}
\alias{slicePNG}
\title{Vectorized NIfTI Slice Rendering with Physical Scaling}
\description{
  Extracts and renders multiple 2D slices from a NIfTI volume into PNG files.
  This function is optimized for speed by unzipping and loading the volume
  once to process all requested sagittal, coronal, and axial slices in a single pass.
  It supports physical scaling (mm-to-pixel), dual-mode thresholding, and the
  generation of reusable annotation templates.
}
\usage{
slicePNG(nii_data, nii_vol = 1, nii_mask = NULL, mask_vol = 1,
         slice_x = NULL, slice_y = NULL, slice_z = NULL,
         bg_color = "black", color = c("black", "white"),
         threshold_pct = NULL, threshold_value = NULL, scale = 1,
         draw_scale = FALSE, draw_side = FALSE, draw_coords = FALSE,
         draw_cbar = NULL, draw_mask = FALSE, draw_label_layer = FALSE,
         file_name = NULL, dir_scratch, dir_save)
}
\arguments{
  \item{nii_data}{Character. Path to the input NIfTI file (.nii or .nii.gz).}
  \item{nii_vol}{Integer. Volume index for 4D NIfTI data (defaults to 1).}
  \item{nii_mask}{Character. Optional path to a NIfTI mask file.}
  \item{mask_vol}{Integer. Volume index for 4D NIfTI masks.}
  \item{slice_x, slice_y, slice_z}{Numeric vectors. Indices for Sagittal (X), Coronal (Y), and Axial (Z) slices to be rendered.}
  \item{bg_color}{Character. Background canvas color. Use \code{"none"} for transparent foreground layers.}
  \item{color}{Character vector. Colors for the intensity color ramp.}
  \item{threshold_pct}{Numeric vector. Percentiles (0-1) for robust scaling.}
  \item{threshold_value}{Numeric vector. Raw intensity values for thresholding. Supports \code{NA} for data-driven bounds.}
  \item{scale}{Numeric. Pixels per millimeter (mm) for physical scaling.}
  \item{draw_scale}{Logical. If TRUE, generates a plane-specific scale bar template (\code{label_[plane]_scale.png}).}
  \item{draw_side}{Logical. If TRUE, generates a plane-specific laterality template (\code{label_[plane]_side.png}) with a directional arrow.}
  \item{draw_coords}{Logical. If TRUE, bakes world-space (mm) coordinates directly onto each slice image.}
  \item{draw_cbar}{Character. Orientation of the color bar: \code{"vertical"} or \code{"horizontal"}.}
  \item{draw_mask}{Logical. If TRUE, generates binary alpha masks (\code{_mask.png}) for each slice.}
  \item{draw_label_layer}{Logical. If TRUE, saves the annotation templates to \code{dir_save} in addition to \code{dir_scratch}.}
  \item{file_name}{Character. Prefix for output files. Defaults to \code{"slice"}.}
  \item{dir_scratch}{Character. Path to local scratch for decompression and temporary files.}
  \item{dir_save}{Character. Directory where final PNGs are saved.}
}
\details{
  \bold{Vectorization:} By providing vectors to \code{slice_x}, \code{slice_y}, or \code{slice_z},
  the function minimizes disk I/O. Background files are unzipped once, and the requested
  volume is loaded into memory once to extract all slices across all three planes.

  \bold{Laterality:} When \code{draw_side = TRUE}, the function calculates the
  anatomical Left/Right orientation using the NIfTI \code{srow_x} header. The resulting
  template includes an arrow (e.g., \eqn{\leftarrow} L) indicating which side of the image
  corresponds to that hemisphere.
}
\value{
  Character vector. File paths to all generated slice PNG images.
}
