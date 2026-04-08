cluster3D <- function(bin_array, connectivity = 26) {
  dims <- dim(bin_array)
  cdim <- cumprod(c(1, dims[1:2]))

  # 1. Linear Offsets (Internal logic remains the same)
  neighborhood <- switch(as.character(connectivity),
                         `6`  = t(matrix(c(1,2,2, 2,1,2, 2,2,1, 3,2,2, 2,3,2, 2,2,3), nrow=3)),
                         `18` = t(matrix(c(1,1,2, 1,2,1, 1,2,2, 1,2,3, 1,3,2, 2,1,1, 2,1,2, 2,1,3, 2,2,1, 2,2,3,
                                           2,3,1, 2,3,2, 2,3,3, 3,1,2, 3,2,1, 3,2,2, 3,2,3, 3,3,2), nrow=3)),
                         `26` = t(matrix(c(1,1,1, 1,1,2, 1,1,3, 1,2,1, 1,2,2, 1,2,3, 1,3,1, 1,3,2, 1,3,3, 2,1,1,
                                           2,1,2, 2,1,3, 2,2,1, 2,2,3, 2,3,1, 2,3,2, 2,3,3, 3,1,1, 3,1,2, 3,1,3,
                                           3,2,1, 3,2,2, 3,2,3, 3,3,1, 3,3,2, 3,3,3), nrow=3)))

  # Calculate offsets once
  center_idx <- sum(c(1,1,1) * c(1, dims[1], dims[1]*dims[2])) # conceptual
  offsets <- apply(neighborhood, 1, function(p) sum((p-1) * c(1, dims[1], dims[1]*dims[2]))) -
    sum((c(2,2,2)-1) * c(1, dims[1], dims[1]*dims[2]))

  # 2. Linearize for memory efficiency
  bin_vec <- as.vector(bin_array)
  connected <- integer(length(bin_vec)) # Uses 4-byte integers instead of 8-byte doubles
  num_clusters <- 0

  # 3. Single Linear Scan
  # Instead of X,Y,Z loops, we just find the indices of 1s
  active_indices <- which(bin_vec == 1)

  for (i in active_indices) {
    if (bin_vec[i] == 1) {
      num_clusters <- num_clusters + 1
      connected[i] <- num_clusters
      bin_vec[i] <- 0

      # Queue for current cluster
      stack <- i
      while (length(stack) > 0) {
        # Pop last element (DFS approach is often more memory stable in R)
        curr <- stack[length(stack)]
        stack <- stack[-length(stack)]

        # Check neighbors
        v_neighbors <- curr + offsets

        # Filter bounds and active voxels
        # This keeps the 'idx' tiny and memory-friendly
        valid <- v_neighbors[v_neighbors > 0 & v_neighbors <= length(bin_vec)]
        hits <- valid[bin_vec[valid] == 1]

        if (length(hits) > 0) {
          connected[hits] <- num_clusters
          bin_vec[hits] <- 0
          stack <- c(stack, hits)
        }
      }
    }
  }

  return(array(connected, dim = dims))
}
