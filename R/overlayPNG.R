overlayPNG <- function(bg_nii,
                       bg_mask,
                       bg_vol = 1,
                       bg_mask_vol = 1,
                       bg_threshold_pct = c(0.025, 0.975),
                       bg_threshold_value,
                       bg_color = c("black", "white"),
                       canvas_color = "black",
                       fg_nii_list,
                       fg_mask_list,
                       fg_vol_list = 1,
                       fg_mask_vol_list = 1,
                       fg_colors = list(c("cyan", "blue"), c("red", "yellow")),
                       fg_threshold_pct,
                       fg_threshold_value,
                       fg_alphas = 1,
                       roi_nii_list,
                       roi_vol_list = 1,
                       roi_value = "all",
                       roi_color = "hotpink",
                       roi_alpha = 1, 
                       roi_outline = TRUE,
                       slice_x = NULL, slice_y = NULL, slice_z = NULL,
                       scale = 1,
                       draw_side = FALSE,
                       draw_coords = FALSE,
                       draw_scale = FALSE,
                       draw_cbar = NULL,
                       apply_labels = FALSE,
                       file_name = NULL,
                       dir_scratch,
                       cleanup = TRUE,
                       dir_save) {
  
  # Setup Local Files ----------------------------------------------------------
  if (missing(bg_mask)) { bg_mask <- NULL }
  n_fg <- 0
  if (!missing(fg_nii_list)) {
    n_fg <- length(fg_nii_list)
    if (!missing(fg_mask_list) & length(fg_mask_list) == 1) {
      fg_mask_list <- rep(fg_mask_list, i)
    } else {
      local_fg_mask <- rep(NA,n_fg)
    }
  }
  if (!missing(roi_nii_list)) {
    n_rois <- length(roi_nii_list)
    if (length(roi_vol_list) != n_rois) { roi_vol_list <- rep(roi_vol_list[1], n_rois) }
    if (length(roi_color) != n_rois) { roi_color <- rep(roi_color, n_rois) }
    if (length(roi_alpha) != n_rois) { roi_alpha <- rep(roi_alpha, n_rois) } # Add this line
  }
  
  # 3. Resolve Background & Setup Canvas ---------------------------------------
  if (missing(bg_threshold_value)) { bg_threshold_value <- NULL }
  bg_paths <- slicePNG(nii_data = bg_nii, nii_mask = bg_mask,
                       nii_vol = bg_vol, mask_vol = bg_mask_vol,
                       slice_x = slice_x, slice_y = slice_y, slice_z = slice_z,
                       scale = scale, color = bg_color, bg_color = canvas_color,
                       threshold_pct = bg_threshold_pct,
                       threshold_value = bg_threshold_value,
                       draw_scale = draw_scale, draw_side = draw_side,
                       draw_coords = draw_coords, draw_label_layer = TRUE,
                       file_name = "bg",
                       dir_scratch = sprintf("%s/slice_tmp", dir_scratch),
                       dir_save = dir_scratch)
  
  # 4. Foreground Stacking -----------------------------------------------------
  if (n_fg > 0) {
    if (missing(fg_threshold_pct) & missing(fg_threshold_value)) {
      fg_threshold_pct <- rep(c(0,1), n_fg)
      fg_threshold_value <- NULL
    } else if (missing(fg_threshold_pct)) {
      fg_threshold_pct <- NULL
      if (is.list(fg_threshold_value) & length(fg_threshold_value) == 1) {
        fg_threshold_value <- rep(fg_threshold_value, n_fg)
      } else {
        fg_threshold_value <- rep(list(fg_threshold_value[1:2]), n_fg)
      }
    } else if (missing(fg_threshold_value)) {
      fg_threshold_value <- NULL
      if (is.list(fg_threshold_pct) & length(fg_threshold_pct) == 1) {
        fg_threshold_pct <- rep(fg_threshold_pct, n_fg)
      } else {
        fg_threshold_pct <- rep(list(fg_threshold_pct[1:2]), n_fg)
      }
    } else {
      fg_threshold_pct <- NULL
      fg_threshold_value <- NULL
    }
    if (length(fg_vol_list) == 1) { fg_vol_list <- rep(fg_vol_list, n_fg) }
    if (length(fg_mask_vol_list) == 1) { fg_mask_vol_list <- rep(fg_mask_vol_list, n_fg) }
    if (length(fg_alphas) == 1) { fg_alphas <- rep(fg_alphas, n_fg) }
    fg_results <- list()
    for (i in 1:n_fg) {
      fg_results[[i]] <- slicePNG(nii_data = fg_nii_list[[i]],
                                  nii_mask = fg_mask_list[[i]],
                                  nii_vol = fg_vol_list[i],
                                  mask_vol = fg_mask_vol_list[i],
                                  slice_x = slice_x, slice_y = slice_y, slice_z = slice_z,
                                  color = fg_colors[[i]], bg_color=canvas_color,
                                  threshold_pct = fg_threshold_pct[[i]],
                                  threshold_value = fg_threshold_value[[i]],
                                  scale = scale, draw_mask = TRUE, draw_cbar = draw_cbar,
                                  file_name = sprintf("fg%d", i),
                                  dir_scratch = sprintf("%s/slice_tmp", dir_scratch),
                                  dir_save = dir_scratch)
    }
  }
  
  # Composite FG on BG ---------------------------------------------------------
  ## Pre-load ROI volumes to avoid hitting the disk repeatedly
  roi_vols <- list()
  if (!missing(roi_nii_list)) {
    roi_vols <- lapply(seq_along(roi_nii_list), function(i) {
      nifti.io::read.nii.volume(prepNII(roi_nii_list[[i]], "roi", dir_scratch), roi_vol_list[i])
    })
  }
  
  # Recreate the master work list so we know the plane/slice for each path
  work_list <- list()
  if(!is.null(slice_x)) work_list <- c(work_list, lapply(slice_x, function(s) list(p="sagittal", s=s)))
  if(!is.null(slice_y)) work_list <- c(work_list, lapply(slice_y, function(s) list(p="coronal", s=s)))
  if(!is.null(slice_z)) work_list <- c(work_list, lapply(slice_z, function(s) list(p="axial", s=s)))
  
  final_paths <- character()
  
  # --- Persistence Check ---
  # Give the file system a split second to catch up with slicePNG's writes
  Sys.sleep(0.5) 
  # Verify that the expected number of bg_paths actually exists
  # This prevents the "index=0; extent=0" error if a read fails
  actual_bg <- file.exists(bg_paths)
  missing_count <- sum(!file.exists(bg_paths))
  if (missing_count > 0) {
    message(sprintf("Waiting for %d delayed files...", missing_count))
    Sys.sleep(1.5) # Extended wait if files are missing
  }
  
  # Color Bar Prep & Export ---
  master_cbar_ready <- NULL
  cbar_files <- list.files(dir_scratch, pattern = "_cbar\\.png$", full.names = TRUE)
  
  if (!is.null(draw_cbar) && length(cbar_files) > 0) {
    # If not baking into slices, export the bars once for the final montage
    if (!apply_labels) {
      file.copy(cbar_files, file.path(dir_save, basename(cbar_files)), overwrite = TRUE)
    }
    # Read unique bars into memory once to avoid disk I/O in the loop
    # We use unique() in case multiple planes generated identical bars
    master_cbar_ready <- magick::image_read(unique(cbar_files))
    master_cbar_ready <- magick::image_background(master_cbar_ready, canvas_color)
  }
  
  for (s_idx in seq_along(bg_paths)) {
    curr_p <- work_list[[s_idx]]$p
    curr_s <- work_list[[s_idx]]$s
    
    # A. Base Anatomy
    img_stack <- magick::image_read(bg_paths[s_idx])
    img_stack <- magick::image_convert(img_stack, type = "truecoloralpha")
    info <- magick::image_info(img_stack) 
    # B. Add Foreground Layers
    if (n_fg > 0) {
      for (l_idx in seq_along(fg_results)) {
        fg_path <- fg_results[[l_idx]][s_idx]
        mask_path <- gsub("\\.png$", "_mask.png", fg_path)
        fg_img <- magick::image_read(fg_path)
        fg_mask <- magick::image_read(mask_path)
        # 1. Spatial Masking
        fg_trans <- magick::image_composite(fg_img, fg_mask, operator = "CopyOpacity")
        # 2. Global Opacity Adjustment
        # This scales the existing alpha channel by your fg_alphas value (0.0 to 1.0)
        # It ensures the 'holes' stay transparent while the 'solid' bits become see-through
        fg_trans <- magick::image_fx(fg_trans, expression = paste0("a*", fg_alphas[l_idx]), channel = "alpha")      
        img_stack <- magick::image_composite(img_stack, fg_trans, operator = "Over")
      }
    }
    
    # C. Add ROI Layers --------------------------------------------------------
    if (!missing(roi_nii_list)) {
      n_roi_files <- length(roi_vols)
      # Normalize inputs to lists if they aren't already
      if (!is.list(roi_value)) roi_value <- rep(list(roi_value), n_roi_files)
      if (!is.list(roi_color)) roi_color <- rep(list(roi_color), n_roi_files)
      if (!is.list(roi_outline)) roi_outline <- rep(list(roi_outline), n_roi_files)
      # Ensure they match the number of files
      roi_value   <- rep(roi_value,   length.out = n_roi_files)
      roi_color   <- rep(roi_color,   length.out = n_roi_files)
      roi_outline <- rep(roi_outline, length.out = n_roi_files)
      
      for (r_idx in seq_along(roi_vols)) {
        vol_data <- roi_vols[[r_idx]]        
        # Determine labels for THIS specific file
        curr_vals <- roi_value[[r_idx]]
        if (length(curr_vals) == 1 && curr_vals == "all") {
          labels_to_plot <- sort(unique(as.vector(vol_data[vol_data > 0])))
        } else {
          labels_to_plot <- as.numeric(curr_vals)
        }
        # Determine colors for THIS specific file
        curr_colors <- roi_color[[r_idx]]
        curr_outline <- roi_outline[[r_idx]]
        
        # Determine labels first to know HOW MANY colors we need
        if (length(roi_value[[r_idx]]) == 1 && roi_value[[r_idx]] == "all") {
          labels_to_plot <- sort(unique(as.vector(vol_data[vol_data > 0])))
        } else {
          labels_to_plot <- as.numeric(roi_value[[r_idx]])
        }
        # DYNAMIC PALETTE CHECK:
        # If the user passed a function (like viridis or heat.colors)
        if (is.function(curr_colors)) {
          # Execute the function to get exactly one color per label
          curr_colors <- curr_colors(length(labels_to_plot))
        }
        
        for (l_idx in seq_along(labels_to_plot)) {
          curr_lab <- labels_to_plot[l_idx]          
          # Mask out only the current label
          lab_mask <- (vol_data == curr_lab) * 1
          r_slice <- switch(curr_p, 
                            "coronal"  = lab_mask[, curr_s, ],
                            "axial"    = lab_mask[, , curr_s],
                            "sagittal" = lab_mask[curr_s, , ])
          if (sum(r_slice) == 0) next          
          # Determine the pixels to color (Outline vs Solid)
          if (curr_outline) {
            interior <- (r_slice > 0) &
              (cbind(r_slice[,-1], 0) > 0) &
              (cbind(0, r_slice[,-ncol(r_slice)]) > 0) &
              (rbind(r_slice[-1,], 0) > 0) &
              (rbind(0, r_slice[-nrow(r_slice),]) > 0)
            roi_mx <- (r_slice > 0 & !interior)
          } else {
            roi_mx <- (r_slice > 0)
          }
          
          if (sum(roi_mx) > 0) {
            roi_hex <- matrix("transparent", nrow = nrow(r_slice), ncol = ncol(r_slice))            
            # Map color: Modulo ensures we loop back if we have fewer colors than labels
            target_col <- curr_colors[((l_idx - 1) %% length(curr_colors)) + 1]            
            roi_hex[roi_mx] <- target_col
            # Render and composite
            # 1. Read and immediately convert to ensure alpha channel exists
            roi_img <- magick::image_read(roi_hex)
            roi_img <- magick::image_convert(roi_img, type = "truecoloralpha")
            # 2. Re-orient and resize as before
            roi_img <- magick::image_flop(magick::image_rotate(roi_img, 270))
            roi_img <- magick::image_resize(roi_img, 
              geometry = magick::geometry_size_pixels(width = info$width, height = info$height, FALSE))
            # 3. Apply alpha scaling
            # Ensure roi_alphas[r_idx] is passed correctly from the setup
            roi_img <- magick::image_fx(roi_img, expression = paste0(roi_alpha[r_idx], "*a"), channel = "alpha")
            # 4. Composite
            img_stack <- magick::image_composite(img_stack, roi_img, operator = "Over")
          }
        }
      }
    }
    
    # D. Add labels if requested
    side_path <- file.path(dir_scratch, sprintf("label_%s_side.png", curr_p))
    scale_path <- file.path(dir_scratch, sprintf("label_%s_scale.png", curr_p))
    if (apply_labels) {
      curr_p <- work_list[[s_idx]]$p
      if (file.exists(side_path)) {
        side_template <- magick::image_read(side_path)
        img_stack <- magick::image_composite(img_stack, side_template, operator = "Over")
      }
      if (file.exists(scale_path)) {
        scale_template <- magick::image_read(scale_path)
        img_stack <- magick::image_composite(img_stack, scale_template, operator = "Over")
      }
    } else {
      if (file.exists(side_path)) {
        file.copy(side_path, file.path(dir_save, basename(side_path)), overwrite = TRUE)
      }
      if (file.exists(scale_path)) {
        file.copy(scale_path, file.path(dir_save, basename(scale_path)), overwrite = TRUE)
      }
    }
    
    # E. Append Color Bar(s) ---------------------------------------------------
    if (apply_labels && !is.null(master_cbar_ready)) {
      info <- magick::image_info(img_stack)
      cbars_list <- list()
      # Center the pre-loaded bars against the current slice dimensions
      for(i in 1:length(master_cbar_ready)) {
        bar_info <- magick::image_info(master_cbar_ready[i])
        if (draw_cbar == "vertical") {
          cbars_list[[i]] <- magick::image_extent(master_cbar_ready[i],
                                                  geometry = sprintf("%dx%d", bar_info$width, info$height),
                                                  gravity = "Center", color = canvas_color)
        } else {
          cbars_list[[i]] <- magick::image_extent(master_cbar_ready[i],
                                                  geometry = sprintf("%dx%d", info$width, bar_info$height),
                                                  gravity = "Center", color = canvas_color)
        }
      }
      # Assemble and attach
      cbars_padded <- do.call(c, cbars_list)
      is_vert <- (draw_cbar == "vertical")
      cbar_strip <- magick::image_append(cbars_padded, stack = !is_vert) # Stack opposite of orientation
      img_stack <- magick::image_append(c(img_stack, cbar_strip), stack = !is_vert)
    }
    
    # Save overlay -------------------------------------------------------------
    # Determine prefix
    if (is.null(file_name) || file_name == "") {
      prefix <- "overlay"
    } else {
      prefix <- tools::file_path_sans_ext(basename(file_name))
    }
    current_base_name <- sprintf("%s_%s_%03d", prefix, curr_p, curr_s)
    out_path <- file.path(dir_save, paste0(current_base_name, ".png"))
    magick::image_write(img_stack, out_path)
    final_paths <- c(final_paths, out_path)
    rm(img_stack)
    if (s_idx %% 10 == 0) gc() 
  }
  if (cleanup && !is.null(dir_scratch) && dir_scratch != "" && dir_scratch != getwd()) {
    unlink(dir_scratch, recursive = TRUE)
  }
  return(final_paths)
}
