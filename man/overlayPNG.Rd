\name{overlayPNG}
\alias{overlayPNG}
\title{Vectorized Multi-Layer NIfTI Compositing}
\description{
  Composites multiple statistical foreground layers and ROI outlines onto a set
  of anatomical background slices. This function utilizes a vectorized
  architecture to process multiple sagittal, coronal, and axial slices in a
  single pass, significantly reducing NIfTI I/O and processing time.
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
           roi_color = "hotpink", slice_x = NULL, slice_y = NULL,
           slice_z = NULL, scale = 1, draw_side = FALSE,
           draw_coords = FALSE, draw_scale = FALSE, draw_cbar = NULL,
           apply_labels = FALSE, file_name = NULL,
           dir_scratch, dir_save)
}
\arguments{
  \item{bg_nii}{Character. Path to the anatomical background NIfTI.}
  \item{bg_mask}{Character. Optional mask for the background image.}
  \item{bg_vol, bg_mask_vol}{Integer. Volume indices for 4D background files.}
  \item{bg_threshold_pct, bg_threshold_value}{Numeric vectors. Intensity scaling bounds for the background.}
  \item{bg_color}{Character vector. Palette for the background anatomy.}
  \item{canvas_color}{Character. Base color for the image canvas (defaults to \code{"black"}).}
  \item{fg_nii_list, fg_mask_list}{Lists of paths for foreground maps and their respective masks.}
  \item{fg_vol_list, fg_mask_vol_list}{Numeric vectors. Volume indices for each foreground layer.}
  \item{fg_colors}{List of character vectors. Palettes for each foreground layer.}
  \item{fg_threshold_pct, fg_threshold_value}{Lists of numeric vectors. Scaling bounds for foreground maps.}
  \item{fg_alphas}{Numeric vector. Opacity (0-1) for each foreground layer.}
  \item{roi_nii_list, roi_vol_list}{Lists of paths and volume indices for ROI files.}
  \item{roi_value}{Character/Numeric. Value to outline, or \code{"all"}.}
  \item{roi_color}{Character vector. Colors for ROI perimeters.}
  \item{slice_x, slice_y, slice_z}{Numeric vectors. Slices to render across the three anatomical planes.}
  \item{scale}{Numeric. Pixels per mm for physical scaling.}
  \item{draw_side}{Logical. If TRUE, creates laterality templates (\eqn{\leftarrow} L/R).}
  \item{draw_coords}{Logical. If TRUE, bakes mm coordinates onto each slice.}
  \item{draw_scale}{Logical. If TRUE, creates scale bar templates.}
  \item{draw_cbar}{Character. Orientation of the color bar: \code{"vertical"} or \code{"horizontal"}.}
  \item{apply_labels}{Logical. If TRUE, composites laterality, scale, and color bars onto every individual slice.}
  \item{file_name}{Character. Prefix for the output files.}
  \item{dir_scratch}{Character. Path to local scratch for intermediate files.}
  \item{dir_save}{Character. Directory where final PNGs are saved.}
}
\details{
  \bold{Stacking Order:} The function assembles images in the following order:
  (1) Anatomical Background, (2) Foreground Statistical Layers (with Alpha masks),
  (3) ROI Outline perimeters, (4) Orientation and Scale Bar templates (if \code{apply_labels = TRUE}),
  and (5) Color Bars.

  \bold{Efficiency:} By taking vectors for \code{slice_x/y/z}, this function
  internalizes the loop over slices. It calls \code{slicePNG} once per NIfTI file
  to generate all required components, then uses \code{magick} in-memory
  compositing to assemble the layers.
}
\value{
  Character vector. Paths to the generated composite PNG images.
}
