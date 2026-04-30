\name{drawMontage}
\alias{drawMontage}
\title{Publication-Ready Brain Montage Generator}
\description{
  A high-level wrapper that automates the generation of complex brain montages.
  It handles hierarchical slice selection (based on results, ROIs, or anatomy),
  multi-layer compositing for every tile, and grid assembly using a flexible
  layout syntax.
}
\usage{
  drawMontage(bg_nii, bg_mask = NULL, bg_vol = 1, bg_mask_vol = 1,
              bg_threshold_pct = c(0.025, 0.975), bg_threshold_value = NULL,
              bg_color = c("black", "white"), canvas_color = "black",
              layout = "1:x,1:y,1:z", edge_clip = 5,
              fg_nii_list = NULL, fg_mask_list = NULL, fg_vol_list = 1,
              fg_mask_vol_list = 1,
              fg_colors = list(c("cyan", "blue"), c("red", "yellow")),
              fg_threshold_pct = NULL, fg_threshold_value = NULL,
              fg_alphas = 1, roi_nii_list = NULL, roi_vol_list = 1,
              roi_value = "all", roi_color = "hotpink", roi_alphas = 1, roi_outline,
              scale = 1,
              draw_side = FALSE, draw_coords = FALSE, draw_scale = FALSE,
              draw_cbar = NULL, cbar_location = "east",
              file_name = "montage", dir_scratch = NULL,
              cleanup = FALSE, dir_save = getwd())
}
\arguments{
  \item{bg_nii}{Character. Path to the anatomical background NIfTI.}
  \item{bg_mask}{Character. Optional background mask.}
  \item{bg_vol, bg_mask_vol}{Integer. Volume indices for 4D background/mask files.}
  \item{bg_threshold_pct, bg_threshold_value}{Numeric vectors. Scaling bounds for background anatomy.}
  \item{bg_color}{Character vector. Palette for background anatomy.}
  \item{canvas_color}{Character. Color for the grid background and padding.}
  \item{layout}{Character. Grid definition string. Use semicolons (\code{;}) for rows and commas (\code{,}) for column blocks. Format: \code{"n:plane"}. Examples: \code{"5:z"} (one row of 5 axial); \code{"5:x;5:y;5:z"} (3 rows, one per plane).}
  \item{edge_clip}{Numeric. Number of slices to skip at the volume edges during auto-selection.}
  \item{fg_nii_list, fg_mask_list}{Lists of paths for statistical overlays and their masks.}
  \item{fg_vol_list, fg_mask_vol_list}{Numeric vectors. Volume indices for each foreground layer.}
  \item{fg_colors}{List of character vectors. Palettes for each foreground layer.}
  \item{fg_threshold_pct, fg_threshold_value}{Lists of thresholds for foreground scaling.}
  \item{fg_alphas}{Numeric vector. Opacity (0-1) for each foreground layer.}
  \item{roi_nii_list, roi_vol_list}{Lists of paths and volume indices for ROI outlining.}
  \item{roi_value}{Character/Numeric. Specific ROI label to outline.}
  \item{roi_color}{Character vector. Colors for ROI outlines.}
  \item{roi_alpha}{Numeric vector. Opacity (0-1) for each ROI layer.}
  \item{roi_outline}{Logical or list. Whether to render ROIs as outlines (\code{TRUE}) or solid blobs (\code{FALSE}).}
  \item{scale}{Numeric. Pixels per mm for physical scaling.}
  \item{draw_side}{Logical. If TRUE, adds an orientation arrow template to the montage.}
  \item{draw_coords}{Logical. If TRUE, bakes mm coordinates onto every slice.}
  \item{draw_scale}{Logical. If TRUE, adds a scale bar template to the montage.}
  \item{draw_cbar}{Character. Style of color bars: \code{"vertical"} or \code{"horizontal"}.}
  \item{cbar_location}{Character. Compass location for the master color bar strip (e.g., \code{"east"}, \code{"southwest"}).}
  \item{file_name}{Character. Filename for the final PNG montage.}
  \item{dir_scratch}{Character. Path for temporary file processing. If \code{NULL}, a temporary directory is created.}
  \item{cleanup}{Logical. If TRUE, deletes all intermediate PNG files in the scratch directory after completion.}
  \item{dir_save}{Character. Directory to save the final montage.}
}
\details{
  \bold{Hierarchical Selection:} Slices are automatically selected based on the most relevant data available:
    (1) Foreground Masks, (2) ROI extent, (3) Background Mask, or (4) Full Image Volume.

  \bold{Grid Assembly:} Individual tiles are generated via vectorized calls to \code{\link{overlayPNG}}.
  Tiles within a row are midpoint-aligned and padded with \code{canvas_color} to maintain a uniform
  grid regardless of anatomical plane dimensions.

  \bold{Anchor Annotations:} Laterality arrows and scale bars are consolidated and "stamped" onto
  the bottom-right corner of the final grid to minimize visual clutter.
}
\value{
  Character. File path to the finished montage PNG.
}
