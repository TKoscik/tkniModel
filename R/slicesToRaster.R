slicesToRaster <- function(in_array, n_row, n_col, orientation) {
  # dimensions: x, y, z
  slice_count <- 0
  row_list <- list()

  for (i in 1:n_row) {
    col_list <- list()
    for (j in 1:n_col) {
      slice_count <- slice_count + 1

      # Extract slice based on orientation
      slice <- switch(orientation,
                      "coronal"  = in_array[, slice_count, ],
                      "axial"    = in_array[, , slice_count],
                      "sagittal" = in_array[slice_count, , ],
                      stop("Invalid orientation. Use 'coronal', 'axial', or 'sagittal'.")
      )
      col_list[[j]] <- slice
    }
    # Combine columns into a single row-strip
    row_list[[i]] <- do.call(rbind, col_list)
  }

  # Combine row-strips into the final large matrix
  # Note: we use cbind and then flip/rotate as needed by the user's plotting preferences
  out_mx <- do.call(cbind, row_list)

  # Convert matrix to long-format data frame (Base R version of melt)
  df <- expand.grid(x = 1:nrow(out_mx), y = 1:ncol(out_mx))
  df$value <- as.vector(out_mx)

  return(df)
}
