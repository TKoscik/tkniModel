\name{drawOverlay}
\alias{drawOverlay}
\title{Generate Multi-Layered Neuroimaging Montage Plots}
\description{
  Creates high-quality "lightbox" or montage visualizations of neuroimaging data.
  The function stacks multiple statistical overlays onto an anatomical background,
  automatically identifies slices containing data, calculates real-world coordinates
  from NIfTI headers, and includes a locator/index plot and color bars.
}
\usage{
  drawOverlay(anat_nii, over_nii, over_vol, over_color,
              mask_nii = "none", mask_vol = 1, roi_nii = NULL,
              roi_val = NULL, roi_color = "#ff64ff",
              orientation = "coronal", idx_slice = TRUE, cbars = TRUE,
              save_dir = getwd(), file_name = "overlay",
              img_format = "png", img_w = NULL, img_unit = "cm",
              img_dpi = 600, save_plot = TRUE, return_plot = FALSE)
}
\arguments{
  \item{anat_nii}{Character. Path to the anatomical background NIfTI file.}
  \item{over_nii}{Character or list. Path(s) to the statistical overlay NIfTI files.}
  \item{over_vol}{Integer or vector. The volume index/indices to extract from the overlay NIfTI(s).}
  \item{over_color}{Character or list. Color ramp colors for each overlay (e.g., \code{"red"} or \code{c("blue", "white", "red")}).}
  \item{mask_nii}{Character or list. Path(s) to mask files for the overlays. Use \code{"none"} to mask by non-zero voxels.}
  \item{mask_vol}{Integer or vector. Volume index/indices for the mask files.}
  \item{roi_nii}{Character. Optional path to an ROI NIfTI to be drawn as a colored outline.}
  \item{roi_val}{Integer vector. Specific values within the ROI NIfTI to outline.}
  \item{roi_color}{Character. Hex code for the ROI outline color.}
  \item{orientation}{Character. Anatomical plane: \code{"coronal"}, \code{"axial"}, or \code{"sagittal"}.}
  \item{idx_slice}{Logical. If TRUE, includes the side-view locator plot.}
  \item{cbars}{Logical. If TRUE, includes color bars for each overlay.}
  \item{save_dir}{Character. Directory to save the final image.}
  \item{file_name}{Character. Base name for the output file.}
  \item{img_format}{Character. File extension (e.g., \code{"png"}, \code{"pdf"}, \code{"tiff"}).}
  \item{img_w}{Numeric. Width of the output image. Height is calculated automatically.}
  \item{img_unit}{Character. Units for dimensions (\code{"cm"}, \code{"in"}, or \code{"mm"}).}
  \item{img_dpi}{Integer. Resolution of the saved image.}
  \item{save_plot}{Logical. Whether to write the file to disk.}
  \item{return_plot}{Logical. Whether to return the ggplot/grob object to the R session.}
}
\details{
  The function performs automated slice selection by scanning for non-empty voxels in the
  provided overlays or ROIs. It converts 3D volumes into 2D raster montages and uses
  NIfTI \code{sform} headers to provide accurate coordinate labels. Multiple overlays
  are blended additively, allowing for the visualization of overlapping effects.
}
\value{
  Invisibly returns the file path of the saved image, or the plot object if
  \code{return_plot = TRUE}.
}
\examples{
  \dontrun{
    drawOverlay(
      anat_nii   = "T1.nii.gz",
      over_nii   = c("tstat_pos.nii.gz", "tstat_neg.nii.gz"),
      over_vol   = c(1, 1),
      over_color = list(c("black", "red", "yellow"), c("black", "blue", "cyan")),
      orientation = "axial",
      file_name   = "Age_Effect_Results"
    )
  }
}
