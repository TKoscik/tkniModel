#' @title Extract Brain Data by Participant BIDS ID with Mask Flexibility
#' @description Pulls cluster-wise means or specific voxel values from participant NIfTI files
#' with support for multi-volume masks and specific label values.
getNIIData <- function(participant_ids,    # Vector of strings: c("sub-01_ses-1", "sub-02_ses-1")
                       nii_data,           # Path to data dir or "path:_suffix.nii.gz"
                       nii_vol = 1,        # Volume index for 4D participant data
                       mask_nii = NULL,    # Path to labeled cluster mask
                       mask_vol = 1,       # New: Volume index for multi-volume masks
                       mask_value = NULL,  # New: Specific label(s) to extract (e.g., c(10, 22))
                       coords = NULL,      # Nx3 matrix of voxel indices
                       dir_scratch = NULL,
                       cleanup = TRUE,
                       verbose = TRUE) {

  # 1. Setup Scratch Directory ------------------------------------------------
  if (is.null(dir_scratch)) dir_scratch <- file.path(tempdir(), "tkni_extract")
  if (!dir.exists(dir_scratch)) dir.create(dir_scratch, recursive = TRUE)

  # 2. File Identification Logic ----------------------------------------------
  nii_parts <- unlist(strsplit(nii_data, ":"))
  search_dir <- nii_parts[1]
  file_suffix <- if (length(nii_parts) > 1) nii_parts[2] else ""

  if (!dir.exists(search_dir)) stop("Search directory does not exist: ", search_dir)

  all_files <- list.files(search_dir, pattern = "\\.nii(\\.gz)?$",
                          full.names = TRUE, recursive = TRUE)

  matched_files <- sapply(participant_ids, function(id_str) {
    pattern <- paste0(id_str, ".*", file_suffix, "\\.nii(\\.gz)?$")
    matches <- all_files[grepl(pattern, all_files)]
    if (length(matches) == 1) return(matches)
    if (length(matches) > 1) {
      warning(sprintf("Multiple matches for %s. Using first: %s", id_str, basename(matches[1])))
      return(matches[1])
    }
    return(NA)
  })

  res_df <- data.frame(participant_id = participant_ids,
                       nii_file = matched_files,
                       stringsAsFactors = FALSE)

  res_df <- res_df[!is.na(res_df$nii_file), ]
  if (nrow(res_df) == 0) stop("Zero NIfTI files matched. Check participant_ids and suffix.")

  # 3. Data Extraction Logic --------------------------------------------------
  if (!is.null(mask_nii)) {
    # 1. Prepare and read the mask
    m_path <- prepNII(mask_nii, "mask", dir_scratch)
    mask_data <- nifti.io::read.nii.volume(m_path, vol.num = mask_vol)
    # 2. Identify unique labels (excluding background 0)
    if (!is.null(mask_value)) {
      # Use only the user-specified labels
      cluster_ids <- mask_value
    } else {
      # Auto-detect all active labels in the mask
      cluster_ids <- sort(unique(as.vector(mask_data[mask_data > 0])))
    }
    # 3. Pre-calculate the linear indices for each label
    # We name the list elements "cluster_#" to drive the dataframe column names
    mask_indices <- lapply(cluster_ids, function(cid) which(mask_data == cid))
    names(mask_indices) <- paste0("cluster_", cluster_ids)

  } else {
    # If no mask, treat the 'coords' as a single ROI
    # We'll default the name to 'cluster_1' or similar for consistency
    coords_mat <- matrix(coords, ncol = 3)
    mask_indices <- list(cluster_coords = coords_mat)
  }

  # 4. Participant Loop --------------------------------------------------------
  if (verbose) message(sprintf("Extracting data for %d clusters...", length(mask_indices)))
  for (i in 1:nrow(res_df)) {
    loc_f <- prepNII(res_df$nii_file[i], "data", dir_scratch)
    part_vol <- nifti.io::read.nii.volume(loc_f, vol.num = nii_vol)
    # Loop through the list of ROIs we prepared in Step 3
    for (c_name in names(mask_indices)) {
      # part_vol[indices] is incredibly fast in R
      vals <- part_vol[mask_indices[[c_name]]]
      res_df[i, c_name] <- mean(vals, na.rm = TRUE)
    }
    if (verbose && i %% 10 == 0) {
      message(sprintf("Processed %d/%d participants...", i, nrow(res_df)))
    }
  }

  # 5. Cleanup ----------------------------------------------------------------
  if (cleanup) {
    if (verbose) message("Cleaning up temporary scratch files...")

    # Identify files created by THIS function call in the scratch folder
    # We look for files starting with 'data_' or 'mask_' as defined in prepNII
    temp_files <- list.files(dir_scratch,
                             pattern = "^(data_|mask_).*",
                             full.names = TRUE)

    if (length(temp_files) > 0) {
      file.remove(temp_files)
    }

    # If the scratch folder was a temporary one created by this function, remove it
    if (grepl("tkni_extract_", dir_scratch) && length(list.files(dir_scratch)) == 0) {
      unlink(dir_scratch, recursive = TRUE)
    }
  }

  if (verbose) message("Extraction complete.")
  return(res_df)
}
