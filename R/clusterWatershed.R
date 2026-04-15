clusterWatershed <- function(nii_pval, 
                             nii_estimate,
                             nii_test,
                             nii_mask, 
                             effect_name = "effect",
                             effect_volume = 1,
                             mask_volume = 1,
                             peak_thresh = 0.001, 
                             extent_thresh = 0.05, 
                             cluster_size = 25,
                             connectivity = 26,
                             do_pos = TRUE, 
                             do_neg = TRUE,
                             save_clusters = TRUE,
                             save_table = TRUE,
                             save_mask = TRUE,
                             dir_save = getwd(),
                             dir_scratch = NULL) {
# 1. Setup & Directory Prep --------------------------------------------------
  if (is.null(dir_scratch)) dir_scratch <- file.path(tempdir(), "watershed_temp")
  if (!dir.exists(dir_scratch)) dir.create(dir_scratch, recursive = TRUE)
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)

  # Prepare files locally for fast I/O
  m_path <- prepNII(nii_mask, "mask", dir_scratch)
  p_path <- prepNII(nii_pval, "pval", dir_scratch)
  e_path <- prepNII(nii_estimate, "est", dir_scratch)
  t_path <- prepNII(nii_test, "test", dir_scratch)
  
  # Load targeted volumes using specified indices
  mask_vol_data <- nifti.io::read.nii.volume(m_path, vol.num = mask_volume)
  pval_vol     <- nifti.io::read.nii.volume(p_path, vol.num = effect_volume)
  est_vol      <- nifti.io::read.nii.volume(e_path, vol.num = effect_volume)
  test_vol     <- nifti.io::read.nii.volume(t_path, vol.num = effect_volume)
  
  dims <- dim(mask_vol_data)
  
  # Retrieve header info for saving outputs
  pixdim <- nifti.io::info.nii(p_path, field = "pixdim")
  orient <- nifti.io::info.nii(p_path, field = "orientation")

  dirs_to_run <- c()
  if (do_pos) dirs_to_run <- c(dirs_to_run, "pos")
  if (do_neg) dirs_to_run <- c(dirs_to_run, "neg")

  for (curr_dir in dirs_to_run) {
    message(sprintf("Processing %s watershed clusters for %s...", curr_dir, effect_name))
    
    # 2. Define Directional Bounding -------------------------------------------
    dir_mask <- if(curr_dir == "pos") (est_vol > 0) else (est_vol < 0)
    valid_extent <- (pval_vol <= extent_thresh & mask_vol_data > 0 & dir_mask)
    valid_extent[is.na(valid_extent)] <- FALSE # Force NAs to FALSE
    valid_peak   <- (pval_vol <= peak_thresh & mask_vol_data > 0 & dir_mask)
    valid_peak[is.na(valid_peak)] <- FALSE # Force NAs to FALSE
    
    if (sum(valid_peak, na.rm = TRUE) == 0) {
      message(sprintf("No significant %s peaks found. Skipping.", curr_dir))
      next
    }

    # Height map (-log10 p) for the "landscape"
    height_vals <- pval_vol
    height_vals[height_vals == 0] <- min(height_vals[height_vals > 0], na.rm = TRUE)
    heights <- -log10(height_vals) * valid_extent

    # 3. Seed Growth (Flooding Logic from PDF) ---------------------------------
    peak_voxels <- which(valid_peak)
    ordered_vox <- peak_voxels[order(heights[peak_voxels], decreasing = TRUE)]
    
    watershed_map <- array(0, dim = dims)
    n_clusters <- 0
    offsets <- if(connectivity == 26) as.matrix(expand.grid(-1:1, -1:1, -1:1))[-14, ] else
               matrix(c(1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1), ncol=3, byrow=TRUE)

    for (v in ordered_vox) {
      v_coord <- arrayInd(v, dims)
      neighbors_coords <- sweep(offsets, 2, v_coord, "+")
      valid_neigh <- rowSums(neighbors_coords > 0 & sweep(neighbors_coords, 2, dims, "<=")) == 3
      neighbors_lin <- (neighbors_coords[valid_neigh,3]-1)*dims[1]*dims[2] + 
                       (neighbors_coords[valid_neigh,2]-1)*dims[1] + 
                        neighbors_coords[valid_neigh,1]
      
      nearby_labels <- unique(watershed_map[neighbors_lin])
      nearby_labels <- nearby_labels[nearby_labels > 0]
      
      if (length(nearby_labels) == 0) {
        n_clusters <- n_clusters + 1
        watershed_map[v] <- n_clusters
      } else {
        watershed_map[v] <- nearby_labels[1]
      }
    }

    # 4. Filter by Cluster Size -----------------------------------------------
    cluster_counts <- as.data.frame(table(raw_clusters = watershed_map[watershed_map > 0]))
    valid_ids <- as.numeric(as.character(cluster_counts$raw_clusters[cluster_counts$Freq >= cluster_size]))
    
    final_map <- array(0, dim = dims)
    if (length(valid_ids) > 0) {
      for (i in seq_along(valid_ids)) {
        final_map[watershed_map == valid_ids[i]] <- i
      }
    }

    # 5. Table Generation (Peak Identification) -------------------------------
    if (save_table && any(final_map > 0)) {
      cl_counts <- as.data.frame(table(raw_clusters = final_map[final_map > 0]))
      
      cl_table <- data.frame(
        cluster_id     = as.numeric(as.character(cl_counts$raw_clusters)),
        size_voxels    = cl_counts$Freq,
        peak_x = NA, peak_y = NA, peak_z = NA,
        peak_estimate  = NA, peak_test = NA, peak_p = NA
      )
      
      for (j in 1:nrow(cl_table)) {
        # Get voxel coordinates for the current watershed label
        idx <- which(final_map == cl_table$cluster_id[j], arr.ind = TRUE)
        
        # Find peak based on absolute test statistic (t-value)
        cluster_t_vals <- abs(test_vol[idx])
        max_idx_in_cluster <- which.max(cluster_t_vals)
        pk_coord <- idx[max_idx_in_cluster, , drop = FALSE]
        
        # Extract data at peak
        cl_table$peak_x[j]        <- pk_coord[1]
        cl_table$peak_y[j]        <- pk_coord[2]
        cl_table$peak_z[j]        <- pk_coord[3]
        cl_table$peak_p[j]        <- pval_vol[pk_coord]
        cl_table$peak_estimate[j] <- est_vol[pk_coord] # Using est_vol from Section 1
        cl_table$peak_test[j]     <- test_vol[pk_coord]
      }
      
      # Optional: Only keep clusters whose peak actually meets the peak_thresh
      cl_table <- cl_table[cl_table$peak_p <= peak_thresh, ]
      
      # Final sort by size
      if (nrow(cl_table) > 0) {
        cl_table <- cl_table[order(cl_table$size_voxels, decreasing = TRUE), ]
        write.csv(cl_table, file.path(dir_save, paste0(pfx, "_table.csv")), row.names = FALSE)
      } else {
        message(sprintf("No %s watershed clusters met the peak threshold in the final check.", curr_dir))
      }
    }
    
    # 6. Saving Results with Effect-Specific Naming ---------------------------
    # Use effect_name in prefix for organized output
    pfx <- sprintf("effect-%s_dir-%s_p%g_cl%d", effect_name, curr_dir, peak_thresh, cluster_size)
    
    if (save_mask && any(final_map > 0)) {
      nifti.io::init.nii(file.path(dir_save, paste0(pfx, "_mask.nii")), dims=dims, pixdim=pixdim, orient=orient)
      nifti.io::write.nii.volume(file.path(dir_save, paste0(pfx, "_mask.nii")), 1, (final_map > 0)*1)
    }
    
    if (save_clusters && any(final_map > 0)) {
      nifti.io::init.nii(file.path(dir_save, paste0(pfx, "_clusters.nii")), dims=dims, pixdim=pixdim, orient=orient)
      nifti.io::write.nii.volume(file.path(dir_save, paste0(pfx, "_clusters.nii")), 1, final_map)
    }
  }
}
