\name{overlayPNG}
\alias{overlayPNG}
\title{Multi-Layer NIfTI Compositing with Automated Masking}
\description{
  Composites multiple statistical foreground layers and ROI outlines onto an
  anatomical background slice. It uses high-speed \code{magick} operators to
  handle transparency, independent volume selection for 4D files, and
  sign-aware color mapping.
}
\usage{
overlayPNG(bg_nii, bg_mask, bg_vol = 1, bg_mask_vol = 1,
           bg_threshold_pct = c(0.025, 0.975), bg_threshold_value,
           bg_color = c("black", "white"), canvas_color = "black",
           fg_nii_list, fg_mask_list, fg_vol_list = 1,
           fg_mask_vol_list = 1,
           fg_colors = list(c("cyan", "blue"), c("red", "yellow")),
           fg_threshold_pct, fg_threshold_value, fg_alphas = 1,
           roi_nii_list, roi_vol_list = 1, roi_value = "all",
           roi_color = "hotpink", plane = "coronal", slice,
           scale = 1, draw_side = FALSE, draw_coords = FALSE,
           draw_scale = FALSE, draw_cbar = NULL, file_name = NULL,
           dir_scratch, dir_save)
}
\arguments{
  \item{bg_nii}{Character. Path to the anatomical background NIfTI.}
  \item{bg_mask}{Character. Optional mask for the background image.}
  \item{bg_vol}{Integer. Volume index for 4D background files.}
  \item{bg_mask_vol}{Integer. Volume index for 4D background masks.}
  \item{bg_threshold_pct}{Numeric vector. Percentiles for background scaling.}
  \item{bg_threshold_value}{Numeric vector. Raw values for background scaling.}
  \item{bg_color}{Character vector. Palette for the background image.}
  \item{canvas_color}{Character. Color of the base canvas (defaults to \code{"black"}).}
  \item{fg_nii_list}{List of characters. Paths to foreground statistical maps.}
  \item{fg_mask_list}{List of characters. Masks for each foreground layer.}
  \item{fg_vol_list}{Numeric vector. Volume indices for each foreground map.}
  \item{fg_mask_vol_list}{Numeric vector. Volume indices for each foreground mask.}
  \item{fg_colors}{List of character vectors. Color palettes for each layer.}
  \item{fg_threshold_pct}{List of numeric vectors. Percentile thresholds for FG layers.}
  \item{fg_threshold_value}{List of numeric vectors. Raw value thresholds for FG layers.}
  \item{fg_alphas}{Numeric vector. Global opacity (0-1) for each foreground layer.}
  \item{roi_nii_list}{List of characters. Paths to ROI files for outlining.}
  \item{roi_vol_list}{Numeric vector. Volume indices for each ROI file.}
  \item{roi_value}{Character or integer. ROI label value to outline, or \code{"all"}.}
  \item{roi_color}{Character vector. Colors for each ROI outline.}
  \item{plane}{Character. Anatomical plane: \code{"coronal"}, \code{"axial"}, or \code{"sagittal"}.}
  \item{slice}{Integer. The index of the slice to extract.}
  \item{scale}{Numeric. Pixels per millimeter (mm) for physical scaling.}
  \item{draw_side}{Logical. If TRUE, adds orientation labels (L/R or A/P).}
  \item{draw_coords}{Logical. If TRUE, adds world-space (mm) coordinates.}
  \item{draw_scale}{Logical. If TRUE, adds a 5mm scale bar.}
  \item{draw_cbar}{Character. Orientation of generated color bars: \code{"vertical"} or \code{"horizontal"}.}
  \item{file_name}{Character. Custom filename for the resulting PNG.}
  \item{dir_scratch}{Character. Path to local scratch for decompression.}
  \item{dir_save}{Character. Directory where the final PNG is saved.}
}
\details{
  \bold{Layering Architecture:} The function initializes a \code{canvas_color}
  base, renders the anatomy via \code{slicePNG}, and then iteratively
  layers foreground clusters and ROI perimeters. Each foreground layer uses
  a dedicated alpha mask to ensure transparency over the background.

  \bold{Annotations:} To prevent obscuration, labels (L/R, mm, scale bar)
  are generated as a separate transparent layer and composited at the
  very end of the stacking process.
}
\value{
  Character. The file path to the generated composite PNG image.
}
