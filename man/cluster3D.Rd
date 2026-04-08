\name{cluster3D}
\alias{cluster3D}
\title{3D Connected Component Labeling}
\description{
  Identifies and labels connected clusters in a 3D binary array using a seed-fill algorithm.
  Each cluster is assigned a unique integer ID.
}
\usage{
  cluster3D(bin_array, connectivity = 26)
}
\arguments{
  \item{bin_array}{A 3D binary array (values of 0 or 1) representing the thresholded voxels.}
  \item{connectivity}{Integer. The neighborhood definition for connectivity.
    Options are \code{6} (faces), \code{18} (faces and edges), or \code{26} (faces, edges, and corners). Default is \code{26}.}
}
\details{
  The function uses a memory-efficient linear scan and stack-based approach to identify connected components.
  It is optimized for large 3D neuroimaging volumes, avoiding the overhead of triple-nested coordinate loops by using linear indexing.
}
\value{
  A 3D integer array of the same dimensions as \code{bin_array}, where each voxel contains the ID of its cluster (0 for background).
}
\examples{
  # Create a toy 3D array with two separate blocks
  arr <- array(0, dim = c(10, 10, 10))
  arr[2:3, 2:3, 2:3] <- 1
  arr[7:8, 7:8, 7:8] <- 1

  # Label the clusters
  clusters <- cluster3D(arr, connectivity = 26)

  # Table of cluster sizes
  table(clusters[clusters > 0])
}
