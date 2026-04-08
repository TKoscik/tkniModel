\name{slicesToRaster}
\alias{slicesToRaster}
\title{Convert 3D Array Slices to a 2D Raster Data Frame}
\description{
  Reshapes a 3D array into a 2D montage matrix based on specified rows and columns,
  then converts it into a long-format data frame for easy plotting with functions
  like \code{ggplot2::geom_tile}.
}
\usage{
  slicesToRaster(in_array, n_row, n_col, orientation)
}
\arguments{
  \item{in_array}{A 3D numeric array (e.g., a NIfTI volume) to be sliced.}
  \item{n_row}{Integer. The number of rows in the resulting montage grid.}
  \item{n_col}{Integer. The number of columns in the resulting montage grid.}
  \item{orientation}{Character. The anatomical plane to extract. Must be one of
    \code{"axial"}, \code{"coronal"}, or \code{"sagittal"}.}
}
\details{
  The function iterates through the 3D volume and tiles slices into a large 2D matrix.
  It uses efficient list-binding to handle large volumes and returns a data frame
  with \code{x}, \code{y}, and \code{value} columns, compatible with standard R
  plotting libraries.
}
\value{
  A data frame with three columns:
    \itemize{
      \item \code{x}: The horizontal pixel coordinate in the 2D montage.
      \item \code{y}: The vertical pixel coordinate in the 2D montage.
      \item \code{value}: The intensity value of the voxel at that location.
    }
}
\examples{
  # Create a dummy 3D brain volume (100x100x12)
  brain_vol <- array(runif(100*100*12), dim = c(100, 100, 12))

  # Create a 3x4 montage of axial slices
  montage_df <- slicesToRaster(brain_vol, n_row = 3, n_col = 4, orientation = "axial")

  # Example plotting with base R
  # plot(montage_df$x, montage_df$y, col = gray(montage_df$value))
}
