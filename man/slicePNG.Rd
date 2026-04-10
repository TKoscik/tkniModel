\name{slicePNG}
\alias{slicePNG}
\title{Physically Scaled NIfTI Slice Rendering with Advanced Masking}
\description{
  Renders a 2D slice from a NIfTI volume into a PNG file using the \code{magick}
  toolbox. This version supports independent data and mask volume selection for 4D
  files, dual-mode thresholding (percentile or raw value), and a robust
  transparency-aware canvas system for high-quality overlays.
}
\usage{
slicePNG(nii_data, nii_vol = 1, nii_mask = NULL, mask_vol = 1,
         plane = "coronal", slice, bg_color = "black",
         color = c("black", "white"), threshold_pct = NULL,
         threshold_value = NULL, scale = 1, draw_scale = FALSE,
         draw_side = FALSE, draw_coords = FALSE, draw_cbar = NULL,
         draw_mask = FALSE, draw_label_layer = FALSE,
         file_name = NULL, dir_scratch, dir_save)
}
\arguments{
  \item{nii_data}{Character. Path to the input NIfTI file (.nii or .nii.gz).}
  \item{nii_vol}{Integer. Volume index for 4D NIfTI data (defaults to 1).}
  \item{nii_mask}{Character. Optional path to a NIfTI mask file.}
  \item{mask_vol}{Integer. Volume index for 4D NIfTI masks (defaults to 1).}
  \item{plane}{Character. Anatomical plane: \code{"coronal"}, \code{"axial"}, or \code{"sagittal"}.}
  \item{slice}{Integer. The index of the slice to extract.}
  \item{bg_color}{Character. Color of the background canvas. Use \code{"none"} for transparency (essential for overlays).}
  \item{color}{Character vector. Colors for the intensity color ramp.}
  \item{threshold_pct}{Numeric vector. Lower and upper percentiles (0-1) for robust scaling (e.g., \code{c(0.025, 0.975)}).}
  \item{threshold_value}{Numeric vector. Raw intensity values for thresholding. Supports \code{NA} for data-driven bounds (e.g., \code{c(2.5, NA)}).}
  \item{scale}{Numeric. Pixels per millimeter (mm) for physical scaling.}
  \item{draw_scale}{Logical. If TRUE, draws a 5mm-stepped scale bar in the bottom-right.}
  \item{draw_side}{Logical. If TRUE, identifies and prints anatomical labels (L/R or A/P) based on sform.}
  \item{draw_coords}{Logical. If TRUE, prints the world-space (mm) coordinate in the top-left.}
  \item{draw_cbar}{Character. Orientation of the generated color bar: \code{"vertical"} or \code{"horizontal"}.}
  \item{draw_mask}{Logical. If TRUE, generates a separate binary alpha mask (\code{_mask.png}).}
  \item{draw_label_layer}{Logical. If TRUE, saves a separate transparent PNG containing only the annotations (\code{_labels.png}).}
  \item{file_name}{Character. Optional custom filename (strips existing extensions internally).}
  \item{dir_scratch}{Character. Path to local scratch for decompression.}
  \item{dir_save}{Character. Directory where output PNGs are saved.}
}
\details{
  \bold{Thresholding:} Exactly one of \code{threshold_pct} or \code{threshold_value}
  should be used. \code{threshold_value} is preferred for statistical maps where
  specific values (e.g., T-scores) are meaningful. If both are \code{NULL},
  the function defaults to the full intensity range of the data.

  \bold{Compositing:} For foreground overlays, set \code{bg_color = "none"}.
  This ensures that the generated PNG only contains data within the mask,
  making it ready for \code{overlayPNG}.
}
\value{
  Character. The file path to the generated primary PNG image.
}
\examples{
\dontrun{
# Axial T-stat slice with 75\% transparency mask
slicePNG(
  nii_data = "stats.nii.gz",
  nii_vol = 5,
  plane = "axial",
  slice = 87,
  bg_color = "none",
  threshold_value = c(2.5, NA),
  draw_mask = TRUE,
  dir_scratch = "/tmp",
  dir_save = "./plots"
)
}
}
