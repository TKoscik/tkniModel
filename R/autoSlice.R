autoSlice <- function(nii_paths,
                      layout = "1:x,1:y,1:z",
                      edge_clip = 5,
                      dir_scratch = tempdir()) {

  # 1. Parse Layout String ----------------------------------------------------
  # Example: "5:x;5:y;5:z" -> NX=5, NY=5, NZ=5
  nx <- ny <- nz <- 0
  rows <- unlist(strsplit(layout, ";"))
  for (row in rows) {
    cols <- unlist(strsplit(row, ","))
    for (col in cols) {
      parts <- unlist(strsplit(col, ":"))
      num <- as.numeric(parts[1])
      plane <- tolower(parts[2])
      if (plane == "x") nx <- nx + num
      if (plane == "y") ny <- ny + num
      if (plane == "z") nz <- nz + num
    }
  }

  # 2. Find Combined Bounding Box ---------------------------------------------
  # We look for the min/max indices where any of the NIfTIs have non-zero data
  xlim <- ylim <- zlim <- c(Inf, -Inf)

  for (path in nii_paths) {
    # Ensure local unzipped file for nifti.io
    local_nii <- file.path(dir_scratch, paste0("auto_", basename(gsub(".gz$", "", path))))
    if (!file.exists(local_nii)) {
      if (grepl("\\.gz$", path)) R.utils::gunzip(path, destname = local_nii, remove = FALSE)
      else file.copy(path, local_nii)
    }

    # Load first volume
    vol <- nifti.io::read.nii.volume(local_nii, 1)
    # Find all non-zero/non-NA voxel indices
    idx <- which(vol != 0 & !is.na(vol), arr.ind = TRUE)

    if (nrow(idx) > 0) {
      xlim <- c(min(xlim[1], min(idx[, 1])), max(xlim[2], max(idx[, 1])))
      ylim <- c(min(ylim[1], min(idx[, 2])), max(ylim[2], max(idx[, 2])))
      zlim <- c(min(zlim[1], min(idx[, 3])), max(zlim[2], max(idx[, 3])))
    }
  }

  # Fallback to full image if no data found
  if (is.infinite(xlim[1])) {
    img_dims <- nifti.io::info.nii(local_nii, "dims")
    xlim <- c(1, img_dims[1]); ylim <- c(1, img_dims[2]); zlim <- c(1, img_dims[3])
  }

  # 3. Calculate Slice Indices (Ported from Bash Page 11-12) -------------------
  calc_indices <- function(lim, n_req, clip) {
    if (n_req <= 0) return(numeric(0))

    # Apply edge clipping
    start <- lim[1] + clip
    stop  <- lim[2] - clip

    # Ensure we don't clip into non-existence
    if (start >= stop) { start <- lim[1]; stop <- lim[2] }

    if (n_req == 1) {
      return(round((start + stop) / 2))
    } else {
      # Generate sequence and center it (Bash logic)
      step <- (stop - start) / (n_req - 1)
      if (step < 1) step <- 1
      # Adjust start/stop to center the step-sequence in the box
      shift <- (stop - start - (step * (n_req - 1))) / 2
      return(round(seq(from = start + shift, by = step, length.out = n_req)))
    }
  }

  slices <- list(
    x = calc_indices(xlim, nx, edge_clip),
    y = calc_indices(ylim, ny, edge_clip),
    z = calc_indices(zlim, nz, edge_clip)
  )

  return(slices)
}
