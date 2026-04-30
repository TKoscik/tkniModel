\name{overlayPNG}
\alias{overlayPNG}
\title{Vectorized Multi-Layer NIfTI Compositing}
\description{
  Composites multiple statistical foreground layers and ROI outlines onto a set 
  of anatomical background slices. This function utilizes a vectorized 
  architecture to process multiple sagittal, coronal, and axial slices in a 
  single pass. Supports multi-label ROI handling with independent color mapping 
  and outline toggles.
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
           roi_color = "hotpink", roi_alpha = 1, roi_outline = TRUE,
           slice_x = NULL, slice_y = NULL, slice_z = NULL, 
           scale = 1, draw_side = FALSE, draw_coords = FALSE, 
           draw_scale = FALSE, draw_cbar = NULL, 
           apply_labels = FALSE, file_name = NULL, 
           dir_scratch, dir_save)
}
\arguments{
  \item{bg_nii}{Character. Path to the anatomical background NIfTI.}
  \item{bg_mask}{Character. Optional mask for the background image.}
  \item{bg_vol, bg_mask_vol}{Integer. Volume indices for 4D background files.}
  \item{bg_threshold_pct, bg_threshold_value}{Numeric vectors. Scaling bounds for background.}
  \item{bg_color}{Character vector. Palette for the background anatomy.}
  \item{canvas_color}{Character. Base color for the image canvas.}
  \item{fg_nii_list, fg_mask_list}{Lists of paths for foreground maps and masks.}
  \item{fg_vol_list, fg_mask_vol_list}{Numeric vectors. Volume indices for each foreground layer.}
  \item{fg_colors}{List of character vectors. Palettes for each foreground layer.}
  \item{fg_threshold_pct, fg_threshold_value}{Lists of thresholds for foreground maps.}
  \item{fg_alphas}{Numeric vector. Opacity (0-1) for each foreground layer.}
  \item{roi_nii_list, roi_vol_list}{Lists of paths and volume indices for ROI files.}
  \item{roi_value}{List or Vector. Labels to extract from each ROI file. Use \code{"all"} to process every unique non-zero label in a file.}
  \item{roi_color}{List or Vector. Colors for each ROI label. If a list of vectors is provided, colors are mapped to labels within each file.}
  \item{roi_alpha}{Numeric vector. Opacity (0-1) for each ROI layer. Matches the length of \code{roi_nii_list}.}
  \item{roi_outline}{List or Logical. Whether to render ROIs as outlines (\code{TRUE}) or solid blobs (\code{FALSE}). Can be a list matching \code{roi_nii_list}.}
  \item{slice_x, slice_y, slice_z}{Numeric vectors. Slices to render across the three anatomical planes.}
  \item{scale}{Numeric. Pixels per mm for physical scaling.}
  \item{draw_side, draw_coords, draw_scale}{Logical. Annotation toggles for orientation, coordinates, and scale bars.}
  \item{draw_cbar}{Character. Orientation of the color bar: \code{"vertical"} or \code{"horizontal"}.}
  \item{apply_labels}{Logical. If TRUE, bakes annotations and color bars onto every individual slice.}
  \item{file_name}{Character. Prefix for the output files.}
  \item{dir_scratch, dir_save}{Character. Paths for temporary and final file storage.}
}
\details{
  \bold{ROI Handling:} The function now supports advanced atlas rendering. If \code{roi_value} 
  is set to \code{"all"}, the function iterates through every unique integer in the NIfTI 
  volume. Colors in \code{roi_color} are mapped sequentially to these labels. 

  \bold{Stacking Order:} Layers are composited in this sequence: (1) Background, 
  (2) Foreground Stats, (3) ROI Layers, (4) Label Templates, (5) Color Bars. 
  Each ROI label is rendered as a separate \code{magick} object to ensure clean 
  overlapping if multiple regions are present on one slice.
}
\value{
  Character vector. Paths to the generated composite PNG images.
}
