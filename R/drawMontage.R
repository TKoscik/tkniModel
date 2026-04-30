drawMontage <- function(# 1. Structural Inputs
  bg_nii,
  bg_mask = NULL,
  bg_vol = 1,
  bg_mask_vol = 1,

  # 2. Thresholding & Appearance
  bg_threshold_pct = NULL,
  bg_threshold_value = NULL,
  bg_color = c("black", "white"),
  canvas_color = "black",

  layout = "1:x,1:y,1:z",
  edge_clip = 5,

  # 3. Foreground Layers (Lists)
  fg_nii_list = NULL,
  fg_mask_list = NULL,
  fg_vol_list = 1,
  fg_mask_vol_list = 1,
  fg_colors = list(c("cyan", "blue"), c("red", "yellow")),
  fg_threshold_pct = NULL,
  fg_threshold_value = NULL,
  fg_alphas = 1,

  # 4. ROI Layers (Lists)
  roi_nii_list = NULL,
  roi_vol_list = 1,
  roi_value = "all",
  roi_color = "hotpink",
  fg_alpha = 1,
  roi_outline = TRUE,

  # 5. Labels & Scale
  scale = 1,
  draw_side = FALSE,
  draw_coords = FALSE,
  draw_scale = FALSE,
  draw_cbar = NULL,
  cbar_location = "east",

  # 6. File & Directory Management
  file_name = "montage",
  dir_scratch = NULL,
  cleanup = FALSE,
  dir_save = getwd()) {

  # Chunk 1: Setup Directories =================================================
  if (is.null(dir_scratch)) {
    dir_scratch <- file.path(tempdir(), paste0("tkni_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  if (!dir.exists(dir_scratch)) dir.create(dir_scratch, recursive = TRUE)
  if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)

  # Chunk 2: Hierarchical Bounding Box Selection ===============================
  # 1. Define Hierarchy - check FG Masks first, then ROIs, then BG mask, then full image extent
  # 2. Call autoSlice - this returns our vectors: slice_indices$x, slice_indices$y, slice_indices$z
  if (!is.null(fg_mask_list) && any(fg_mask_list != "none")) {
    search_paths <- fg_mask_list[fg_mask_list != "none"]
    message("autoSlice: Using Foreground Masks for bounding box.")

  } else if (!is.null(roi_nii_list)) {
    search_paths <- roi_nii_list
    message("autoSlice: Using ROI files for bounding box.")

  } else if (!is.null(bg_mask)) {
    search_paths <- list(bg_mask)
    message("autoSlice: Using Background Mask for bounding box.")

  } else {
    search_paths <- list(bg_nii)
    message("autoSlice: Using Background Image for bounding box.")
  }
  slice_indices <- autoSlice(nii_paths = search_paths,
                             layout = layout,
                             edge_clip = edge_clip,
                             dir_scratch = dir_scratch)

  # Chunk 3. Generate all Composites via vectorized overlayPNG =================
  # We pass all slice vectors at once to minimize I/O
  all_slice_paths <- overlayPNG(bg_nii = bg_nii,
             bg_mask = bg_mask,
             bg_threshold_pct = bg_threshold_pct,
             bg_threshold_value = bg_threshold_value,
             bg_color = bg_color,
             canvas_color = canvas_color,
             fg_nii_list = fg_nii_list,
             fg_mask_list = fg_mask_list,
             fg_vol_list = fg_vol_list,
             fg_mask_vol_list = fg_mask_vol_list,
             fg_colors = fg_colors,
             fg_threshold_pct = fg_threshold_pct,
             fg_threshold_value = fg_threshold_value,
             fg_alphas = fg_alphas,
             roi_nii_list = roi_nii_list,
             roi_vol_list = roi_vol_list,
             roi_value = roi_value,
             roi_color = roi_color,
             roi_alpha = roi_alpha,
             roi_outline = roi_outline,
             slice_x = slice_indices$x,
             slice_y = slice_indices$y,
             slice_z = slice_indices$z,
             scale = scale,
             draw_side = draw_side,
             draw_coords = draw_coords,
             draw_scale = draw_scale,
             draw_cbar = draw_cbar,
             apply_labels = FALSE,
             dir_scratch = sprintf("%s/overlay", dir_scratch),
             dir_save = dir_scratch)

  # Chunk 4: Grid Assembly & Master Labeling ===================================
 # 1. Split the flat slice_paths into plane-specific pools
  # slice_indices contains the counts of what was generated
  x_count <- length(slice_indices$x)
  y_count <- length(slice_indices$y)
  z_count <- length(slice_indices$z)
  # Map paths to pools (Order in slicePNG: X, then Y, then Z)
  all_imgs <- magick::image_read(all_slice_paths)
  pools <- list(
    x = if(x_count > 0) all_imgs[1:x_count] else NULL,
    y = if(y_count > 0) all_imgs[(x_count + 1):(x_count + y_count)] else NULL,
    z = if(z_count > 0) all_imgs[(x_count + y_count + 1):(x_count + y_count + z_count)] else NULL
  )

  # 2. Tracking pointers for each plane pool
  pool_ptrs <- list(x = 1, y = 1, z = 1)
  row_defs <- strsplit(unlist(strsplit(layout, ";")), ",")
  row_imgs <- list()
  for (i in seq_along(row_defs)) {
    col_imgs <- list()
    for (j in seq_along(row_defs[[i]])) {
      # Parse "5:x"
      def <- unlist(strsplit(row_defs[[i]][j], ":"))
      n_slices <- as.numeric(def[1])
      plane <- def[2]
      # Grab slices from the correct POOL
      start <- pool_ptrs[[plane]]
      end   <- start + n_slices - 1
      block <- pools[[plane]][start:end]
      # Update only the pointer for that specific plane
      pool_ptrs[[plane]] <- end + 1
      # --- Corner Labeling (L/R and Scale) ---
      # [Apply side label to block[[1]] if bottom row and j == 1]
      # [Apply scale bar to block[[n_slices]] if bottom row and j == last]
      if (i == length(row_defs) & j == length(row_defs[[i]])) {
        # 1. Apply Top-Right Laterality Label (e.g., "<- L")
        side_template_path <- list.files(dir_scratch, pattern = "label_.*_side.png", full.names = TRUE)[1]
        if (!is.na(side_template_path)) {
          side_template <- magick::image_read(side_template_path)
          # Resize template to fit current slice dimensions
          side_template <- magick::image_resize(side_template,
             magick::geometry_size_pixels(magick::image_info(block[n_slices])$width,
                                          magick::image_info(block[n_slices])$height, FALSE))
          block[n_slices] <- magick::image_composite(block[n_slices], side_template, operator = "Over")
        }

        # 2. Apply Bottom-Right Scale Bar
        scale_template_path <- list.files(dir_scratch, pattern = "label_.*_scale.png", full.names = TRUE)[1]
        if (!is.na(scale_template_path)) {
          scale_template <- magick::image_read(scale_template_path)
          scale_template <- magick::image_resize(scale_template,
            magick::geometry_size_pixels(magick::image_info(block[n_slices])$width,
                                         magick::image_info(block[n_slices])$height, FALSE))
          block[n_slices] <- magick::image_composite(block[n_slices], scale_template, operator = "Over")
        }
      }
      col_imgs[[j]] <- lapply(seq_along(block), function(idx) block[idx])
    }
    # Flatten the list of lists into a single list of all images for this row
    row_images <- unlist(col_imgs, recursive = FALSE)
    row_images_seq <- do.call(c, row_images)

    # 1. FIND TALLEST IMAGE in this row
    max_h <- max(magick::image_info(row_images_seq)$height)

    # 2. ALIGN AT MIDPOINT and PAD
    # Every image in the row is now forced to max_h, centered, with canvas_color
    centered_row_list <- lapply(row_images, function(img) {
      curr_info <- magick::image_info(img)
      magick::image_extent(img,
                           geometry = sprintf("%dx%d", curr_info$width, max_h),
                           gravity = "Center",
                           color = canvas_color)
    })

    # 3. Create the row strip from centered images
    centered_row_seq <- do.call(c, centered_row_list)
    row_strip <- magick::image_append(centered_row_seq, stack = FALSE)

    # Fill any remaining row-level gaps
    row_imgs[[i]] <- magick::image_background(row_strip, canvas_color)
  }

  # Final vertical stack
  final_montage <- magick::image_append(do.call(c, row_imgs), stack = TRUE)
  final_montage <- magick::image_background(final_montage, canvas_color)
  # Chunk 5: Master Color Bar Assembly =========================================
  if (!is.null(draw_cbar)) {
    all_cbars <- list.files(dir_scratch, pattern = "_cbar\\.png$", full.names = TRUE)
    if (length(all_cbars) > 0) {
      cbar_imgs <- magick::image_read(unique(all_cbars))
      cbar_imgs <- magick::image_background(cbar_imgs, canvas_color)

      # Logic:
      # If style is "horizontal", stack cbars top-to-bottom (stack = TRUE)
      # If style is "vertical", stack cbars left-to-right (stack = FALSE)
      is_stacking_vertically <- (draw_cbar == "horizontal")
      master_cbar <- magick::image_append(cbar_imgs, stack = is_stacking_vertically)

      # 2. Add "Room" to the montage (Padding)
      b_info <- magick::image_info(master_cbar)
      m_info <- magick::image_info(final_montage)
      loc <- tolower(cbar_location)

      # Expand canvas based on location to prevent covering data
      if (grepl("south", loc)) {
        final_montage <- magick::image_extent(final_montage,
                                              geometry = sprintf("%dx%d", m_info$width, m_info$height + b_info$height),
                                              gravity = "North", color = canvas_color)
      } else if (grepl("north", loc)) {
        final_montage <- magick::image_extent(final_montage,
                                              geometry = sprintf("%dx%d", m_info$width, m_info$height + b_info$height),
                                              gravity = "South", color = canvas_color)
      } else if (grepl("east", loc)) {
        final_montage <- magick::image_extent(final_montage,
                                              geometry = sprintf("%dx%d", m_info$width + b_info$width, m_info$height),
                                              gravity = "West", color = canvas_color)
      } else if (grepl("west", loc)) {
        final_montage <- magick::image_extent(final_montage,
                                              geometry = sprintf("%dx%d", m_info$width + b_info$width, m_info$height),
                                              gravity = "East", color = canvas_color)
      }

      # 3. Composite without scaling
      final_montage <- magick::image_composite(final_montage, master_cbar,
                                               gravity = cbar_location,
                                               operator = "Over")
    }
  }

  # Chunk 6: Save and Finalize =================================================
  out_path <- file.path(dir_save, paste0(file_name, ".png"))
  magick::image_write(final_montage, out_path)

  # Chunk 7: Cleanup ===========================================================
  # We delete only the PNGs we created in the scratch directory
  # to keep the folder clean for the next run.
  temp_files <- list.files(dir_scratch, pattern = "\\.png$", full.names = TRUE)
  if (length(temp_files) > 0) { file.remove(temp_files) }
  if (cleanup) { unlink(dir_scratch, recursive = TRUE) }
  message(sprintf("Montage successfully saved to: %s", out_path))
  return(out_path)
}
