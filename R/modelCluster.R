modelCluster <- function(nii_estimate,
                         nii_test,
                         nii_pval,
                         nii_mask        = NULL,
                         effect_name     = "effect",
                         effect_volume   = 1,
                         mask_volume     = 1,
                         roi_thresh      = 0.05,
                         peak_thresh     = 0.001,
                         cluster_size    = 25,
                         connectivity    = 26L,
                         save_clusters   = TRUE,
                         save_table      = TRUE,
                         save_mask       = TRUE,
                         do_pos          = TRUE,
                         do_neg          = TRUE,
                         dir_save        = getwd(),
                         dir_scratch     = NULL,
                         verbose         = TRUE,
                         cleanup         = TRUE) {

  # 1. Validation -------------------------------------------------------------
  missing_input <- FALSE
  if (missing(nii_test)) { message("Error: 'nii_test' is required."); missing_input <- TRUE }
  if (missing(nii_pval)) { message("Error: 'nii_pval' is required."); missing_input <- TRUE }
  if (missing_input) stop("Missing required arguments, aborting.", call. = FALSE)

  # 2. Setup Directories ------------------------------------------------------
  if (is.null(dir_scratch)) {
    timestamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
    dir_scratch <- file.path(tempdir(), paste0("tkni_clusterEffect_", timestamp))
  }

  if (!dir.exists(dir_scratch)) dir.create(dir_scratch, recursive = TRUE)
  if (!dir.exists(dir_save))    dir.create(dir_save, recursive = TRUE)

  # 3. Load Libraries ---------------------------------------------------------
  libs <- c("nifti.io", "R.utils", "tools")
  for (lib in libs) {
    if (!require(lib, character.only = TRUE, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", lib))
    }
  }

  # 4. Copy and Decompress files to scratch -----------------------------------
  # Helper to handle transfer and decompression
  prep_nii <- function(orig_path, target_name, scratch_path) {
    if (is.null(orig_path)) return(NULL)
    dest <- file.path(scratch_path, paste0(target_name, ".nii"))
    if (tools::file_ext(orig_path) == "gz") {
      R.utils::gunzip(filename = orig_path, destname = dest, remove = FALSE, overwrite = TRUE)
    } else {
      file.copy(from = orig_path, to = dest, overwrite = TRUE)
    }
    return(dest)
  }

  nii_estimate <- prep_nii(nii_estimate, "estimate", dir_scratch)
  nii_test     <- prep_nii(nii_test,     "test",     dir_scratch)
  nii_pval     <- prep_nii(nii_pval,     "pval",     dir_scratch)
  nii_mask     <- prep_nii(nii_mask,     "mask",     dir_scratch)

  # 5. Load Image Metadata ----------------------------------------------------
  # Using nii_estimate as the reference for headers
  img_dims <- nifti.io::info.nii(nii_test, "dims")
  pixdim   <- nifti.io::info.nii(nii_test, "pixdim")
  orient   <- nifti.io::info.nii(nii_test, "orient")

  # 6. Load and Process Volumes -----------------------------------------------
  # Load mask (if none provided, create a full-brain mask of 1s)
  if (!is.null(nii_mask)) {
    mask <- nifti.io::read.nii.volume(nii_mask, mask_volume)
  } else {
    mask <- array(1, dim = img_dims[1:3])
  }
  mask <- (mask > 0) * 1 # Binarize

  # Load Estimate (Beta/Coefficient)
  estimate_vol <- nifti.io::read.nii.volume(nii_estimate, effect_volume)
  estimate_vol[is.na(estimate_vol)] <- 0
  estimate_vol <- estimate_vol * mask

  # Load Test Statistic (T/F map)
  test_vol <- nifti.io::read.nii.volume(nii_test, effect_volume)
  test_vol[is.na(test_vol)] <- 0
  test_vol <- test_vol * mask

  # Load P-values
  pval_vol <- nifti.io::read.nii.volume(nii_pval, effect_volume)
  pval_vol[is.na(pval_vol)] <- 0
  pval_vol <- pval_vol * mask

  # 7. Processing Loop (Positive and Negative Effects) ------------------------
  for (side in c("pos", "neg")) {
    # Skip if the user disabled this direction
    if (side == "pos" && !do_pos) next
    if (side == "neg" && !do_neg) next

    if (verbose) message(sprintf("Processing %s clusters...", side))

    # 1. Mask by sign and threshold by p-value
    if (side == "pos") {
      sig_mask <- (pval_vol < roi_thresh) * (estimate_vol > 0) * 1
    } else {
      sig_mask <- (pval_vol < roi_thresh) * (estimate_vol < 0) * 1
    }

    # 2. Generate clusters using our cluster3D function
    raw_clusters <- cluster3D(sig_mask, connectivity = connectivity)

    # 3 & 4. Initial cluster table and sorting
    cl_counts <- as.data.frame(table(raw_clusters))
    cl_counts <- cl_counts[cl_counts$raw_clusters != 0, ] # Remove background

    if (nrow(cl_counts) == 0) {
      if (verbose) message(sprintf("No %s clusters found at p < %g.", side, roi_thresh))
      next
    }

    # 5. Threshold by cluster size
    cl_counts$Freq <- as.numeric(cl_counts$Freq)
    cl_counts <- cl_counts[cl_counts$Freq >= cluster_size, ]

    if (nrow(cl_counts) == 0) {
      if (verbose) message(sprintf("No %s clusters survived size threshold (%d).", side, cluster_size))
      next
    }

    # 6. Calculate peak values and coordinates
    # Initialize table with your required columns
    cl_table <- data.frame(
      cluster_id     = as.numeric(as.character(cl_counts$raw_clusters)),
      size_voxels    = cl_counts$Freq,
      peak_x = NA, peak_y = NA, peak_z = NA,
      peak_estimate  = NA, peak_test = NA, peak_p = NA
    )

    for (j in 1:nrow(cl_table)) {
      idx <- which(raw_clusters == cl_table$cluster_id[j], arr.ind = TRUE)

      # Find peak based on absolute test statistic
      # Using local indices for the test_vol subset
      cluster_vals <- abs(test_vol[idx])
      max_idx_in_cluster <- which.max(cluster_vals)
      pk_coord <- idx[max_idx_in_cluster, , drop = FALSE]

      # Extract values at peak coordinate
      cl_table$peak_x[j]        <- pk_coord[1]
      cl_table$peak_y[j]        <- pk_coord[2]
      cl_table$peak_z[j]        <- pk_coord[3]
      cl_table$peak_p[j]        <- pval_vol[pk_coord]
      cl_table$peak_estimate[j] <- estimate_vol[pk_coord]
      cl_table$peak_test[j]     <- test_vol[pk_coord]
    }

    # 7. Threshold by peak p-value
    cl_table <- cl_table[cl_table$peak_p < peak_thresh, ]

    if (nrow(cl_table) == 0) {
      if (verbose) message(sprintf("No %s clusters survived peak p threshold (%g).", side, peak_thresh))
      next
    }

    # 8. Reorder by size and renumber cluster map
    cl_table <- cl_table[order(-cl_table$size_voxels), ]
    final_cluster_map <- array(0, dim = img_dims[1:3])

    for (new_id in 1:nrow(cl_table)) {
      vox_idx <- which(raw_clusters == cl_table$cluster_id[new_id])
      final_cluster_map[vox_idx] <- new_id
    }

    # Update table IDs to match the new map
    cl_table$cluster_id <- 1:nrow(cl_table)

    # 9. Output results
    tpfx <- file.path(dir_save, sprintf("effect-%s_dir-%s_p-%0.0g_peak-%0.0g_cl-%d",
                                        effect_name, side, roi_thresh, peak_thresh, cluster_size))

    if (save_clusters) {
      out_nii <- paste0(tpfx, "_cluster.nii")
      nifti.io::init.nii(out_nii, dims = img_dims, pixdim = pixdim, orient = orient)
      nifti.io::write.nii.volume(out_nii, 1, final_cluster_map)
      R.utils::gzip(out_nii, overwrite = TRUE)
    }

    if (save_mask) {
      out_mask <- paste0(tpfx, "_mask.nii")
      nifti.io::init.nii(out_mask, dims = img_dims, pixdim = pixdim, orient = orient)
      nifti.io::write.nii.volume(out_mask, 1, (final_cluster_map > 0) * 1)
      R.utils::gzip(out_mask, overwrite = TRUE)
    }

    if (save_table) {
      write.csv(cl_table, paste0(tpfx, "_table.csv"), row.names = FALSE)
    }
  }

  # Final Cleanup
  if (cleanup) {
    unlink(dir_scratch, recursive = TRUE)
    if (verbose) message("Scratch directory cleared.")
  }

  return(invisible(dir_save))
}
