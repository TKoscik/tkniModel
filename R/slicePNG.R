slicePNG <- function(nii_data, nii_vol = 1,
                     nii_mask = NULL, mask_vol = 1,
                     slice_x = NULL, slice_y = NULL, slice_z = NULL,
                     bg_color = "black",
                     color = c("black", "white"),
                     threshold_pct = NULL,   # e.g., c(0.025, 0.975)
                     threshold_value = NULL, # e.g., c(2.5, 5.0)
                     scale = 1,
                     draw_scale = FALSE,
                     draw_side = FALSE,
                     draw_coords = FALSE,
                     draw_cbar = NULL,
                     draw_mask = FALSE,
                     draw_label_layer = FALSE,
                     file_name = NULL,
                     dir_scratch,
                     dir_save,
                     cleanup = TRUE) {

  if (!is.null(threshold_pct) && !is.null(threshold_value)) {
    stop("Provide either threshold_pct OR threshold_value, not both.")
  }
  if (is.null(slice_x) & is.null(slice_y) & is.null(slice_z)) {
    stop("You must provide slices to plot.")
  }

  # 1. Setup Scratch and Local Files ------------------------------------------
  if (!dir.exists(dir_scratch)) dir.create(dir_scratch, recursive = TRUE)
  local_nii <- prepNII(nii_data, "data", dir_scratch)
  local_mask <- prepNII(nii_mask, "mask", dir_scratch)

  # 2. Get Metadata and Thresholds --------------------------------------------
  spacing <- nifti.io::info.nii(local_nii, "spacing")
  img_vol <- nifti.io::read.nii.volume(local_nii, nii_vol)
  dims <- dim(img_vol)

  # 3. Load mask ---------------------------------------------------------------
  if (!is.null(local_mask)) {
    m_vol <- nifti.io::read.nii.volume(local_mask, mask_vol)
    m_vol <- ((m_vol != 0) & (!is.na(img_vol)))*1
  } else {
    m_vol <- img_vol * 0 + 1
  }

  # 4. Resolve Thresholds ------------------------------------------------------
  if(!is.null(local_mask)) {
    vals <- img_vol[m_vol != 0]
  } else {
    vals <- img_vol[!is.na(img_vol)]
  }
  v_min <- min(vals, na.rm=TRUE)
  v_max <- max(vals, na.rm=TRUE)
  if (!is.null(threshold_pct)) {
    v_bounds <- quantile(vals, probs = threshold_pct, na.rm = TRUE)
    v_min <- v_bounds[1]
    v_max <- v_bounds[2]
  } else if (!is.null(threshold_value)) {
    if (!is.na(threshold_value[1])) { v_min <- threshold_value[1] }
    if (!is.na(threshold_value[2])) { v_max <- threshold_value[2] }
  }

  if (is.null(local_mask)) {
    m_vol[img_vol < min(c(v_min, v_max)) & img_vol > max(c(v_min, v_max))] <- 0
  }
  
  # 5. Create the slice extraction loop ----------------------------------------
  work_list <- list()
  if(!is.null(slice_x)) work_list <- c(work_list, lapply(slice_x, function(s) list(p="sagittal", s=s)))
  if(!is.null(slice_y)) work_list <- c(work_list, lapply(slice_y, function(s) list(p="coronal", s=s)))
  if(!is.null(slice_z)) work_list <- c(work_list, lapply(slice_z, function(s) list(p="axial", s=s)))
  plane_to_arg <- list(sagittal = "slice_x", coronal = "slice_y", axial = "slice_z")
  output_paths <- character()

  # calculate dimensions for each plane -----------------------------------------
  p_dims <- list(
    sagittal = list(w = dims[2] * spacing[2], h = dims[3] * spacing[3]),
    coronal  = list(w = dims[1] * spacing[1], h = dims[3] * spacing[3]),
    axial    = list(w = dims[1] * spacing[1], h = dims[2] * spacing[2])
  )
  # Calculate font sizing based on scale
  base_font <- max(10, round(10 * scale))
  padding <- 15

  # 5. Generate Reusable Templates (Side and Scale) ----------------------------
  ## Scale Labels
  if (draw_scale) {
    for (p in names(p_dims)) {
      arg_name <- plane_to_arg[[p]]
      if (!is.null(get(arg_name))) {
        tw <- round(p_dims[[p]]$w * scale); th <- round(p_dims[[p]]$h * scale)
        # Calculate bar physics
        target_mm <- p_dims[[p]]$w * 0.15
        scale_len_mm <- round(target_mm / 5) * 5
        if (scale_len_mm == 0) scale_len_mm <- 5
        bar_w_px <- scale_len_mm * scale
        bar_h_px <- max(1, round(0.3 * scale))
        c_img <- magick::image_blank(tw, th, color = "none")
        c_img <- magick::image_draw(c_img) # This is the missing step!
        symbols(tw - padding - (bar_w_px/2),
                th - padding - (base_font * 1.5) - (bar_h_px/2),
                rectangles = matrix(c(bar_w_px, bar_h_px), ncol=2),
                inches = FALSE, add = TRUE, bg = "white", fg = NA)
        dev.off()
        c_img <- magick::image_annotate(c_img, sprintf("%g mm", scale_len_mm),
                                        gravity = "southeast", color = "white", size = base_font,
                                        location = sprintf("+%d+%d", round(padding + (bar_w_px/2) - (base_font)), padding))
        magick::image_write(c_img, file.path(dir_scratch, sprintf("label_%s_scale.png", p)))
        list.files(dir_scratch)
        if (draw_label_layer) {
          magick::image_write(c_img, file.path(dir_save, sprintf("label_%s_scale.png", p)))
        }
      }
    }
  }

  ## Side labels
  if (draw_side) {
    srow_x <- unlist(nifti.io::info.nii(local_nii, "srow_x"))
    side_char <- if (srow_x[1] > 0) "R" else "L"
    label_text <- paste0("\u2190 ", side_char) # "<- L" or "<- R"
    for (p in c("coronal", "axial")) {
      arg_name <- plane_to_arg[[p]]
      if (!is.null(get(arg_name))) {
        tw <- round(p_dims[[p]]$w * scale)
        th <- round(p_dims[[p]]$h * scale)
        side_img <- magick::image_blank(tw, th, color = "none")
        side_img <- magick::image_annotate(side_img, label_text,
          gravity = "northeast", color = "white", size = base_font * 1.5,
          location = sprintf("+%d+%d", padding, padding))
        magick::image_write(side_img, file.path(dir_scratch, sprintf("label_%s_side.png", p)))
        if (draw_label_layer) {
          magick::image_write(side_img, file.path(dir_save, sprintf("label_%s_side.png", p)))
        }
      }
    }
  }

  # Loop over slice generation --------------------------------------------------
  for (item in work_list) {
    curr_plane <- item$p
    curr_slice <- item$s

    # A. Extract Slice and Calculate Physical Dimensions ------------------------
    if (curr_plane == "coronal") {
      slice_data <- img_vol[, curr_slice, ]
      phys_w <- dims[1] * spacing[1]; phys_h <- dims[3] * spacing[3]
    } else if (curr_plane == "axial") {
      slice_data <- img_vol[, , curr_slice]
      phys_w <- dims[1] * spacing[1]; phys_h <- dims[2] * spacing[2]
    } else {
      slice_data <- img_vol[curr_slice, , ]
      phys_w <- dims[2] * spacing[2]; phys_h <- dims[3] * spacing[3]
    }

    # B. Generate Main Image ---------------------------------------------------
    ## This is your RGB hex matrix (e.g., "#FF0000")
    clamped_data <- slice_data
    clamped_data[clamped_data < v_min] <- v_min
    clamped_data[clamped_data > v_max] <- v_max
    norm_data <- (clamped_data - v_min) / (v_max - v_min)
    pal <- grDevices::colorRampPalette(color)(256)
    color_indices <- findInterval(norm_data, seq(0, 1, length.out = 256), all.inside = TRUE)
    hex_mx <- matrix(pal[color_indices], nrow = nrow(slice_data))
    hex_mx[is.na(hex_mx)] <- "#000000"

    # C. Build the Alpha Channel from the Mask ---------------------------------
    ## Create alpha values: "FF" (opaque) for data, "00" (transparent) for background
    m_slice <- switch(curr_plane,
                      "sagittal" = m_vol[curr_slice, , ],
                      "coronal"  = m_vol[, curr_slice, ],
                      "axial"    = m_vol[, , curr_slice])
    alpha_mx <- matrix("00", nrow = nrow(slice_data), ncol = ncol(slice_data))
    alpha_mx[m_slice != 0] <- "FF"

    # D. Combine RGB + Alpha into a single RGBA Matrix -------------------------
    # This results in hex codes like "#FF0000FF" (Opaque Red) or "#FF000000" (Transparent)
    rgba_mx <- matrix(paste0(hex_mx, alpha_mx), nrow = nrow(slice_data))
    # Create Magick Image and Composite over bg_color, brain data with transparency ALREADY BUILT IN
    data_img <- magick::image_read(rgba_mx)
    data_img <- magick::image_flop(magick::image_rotate(data_img, 270))

    # E. Resize Image ----------------------------------------------------------
    target_w <- round(phys_w * scale)
    target_h <- round(phys_h * scale)
    data_img <- magick::image_resize(data_img,
      geometry = magick::geometry_size_pixels(width = target_w,
                                              height = target_h,
                                              preserve_aspect = FALSE),
                                    filter = "Catrom")

    # F. Create the final canvas and layer the pre-masked data over it ---------
    bg_canvas <- magick::image_blank(target_w, target_h, color = bg_color)
    img <- magick::image_composite(bg_canvas, data_img, operator = "Over")

    # G. Filenames -------------------------------------------------------------
    if (is.null(file_name)) {
      base_name <- sprintf("slice_%s_%03d", curr_plane, curr_slice)
    } else {
      base_name <- gsub("\\.png$", "", file_name)
    }

    # I. Add Annotations -------------------------------------------------------
    if (!draw_label_layer) {
      ## add scale annotation
      if (draw_scale) {
        scale_tmp <- magick::image_read(file.path(dir_scratch, sprintf("label_%s_scale.png", curr_plane)))
        img <- magick::image_composite(img, scale_tmp, operator = "Over")
      }
      ## add side annotation
      if (draw_side & curr_plane %in% c("coronal", "axial")) {
        side_tmp <- magick::image_read(file.path(dir_save, sprintf("label_%s_scale.png", curr_plane)))
        img <- magick::image_composite(img, side_tmp, operator = "Over")
      }
    }
    ## Add slice Location (World mm)
    if (draw_coords) {
      hdr_field <- switch(curr_plane, "sagittal"="srow_x", "coronal"="srow_y", "axial"="srow_z")
      step_idx <- switch(curr_plane, "sagittal"=1, "coronal"=2, "axial"=3)
      tf <- unlist(nifti.io::info.nii(local_nii, hdr_field))
      world_mm <- round((curr_slice - 1) * tf[step_idx] + tf[4], 1)
      img <- magick::image_annotate(img, sprintf("%g mm", world_mm),
        gravity = "northwest", location = sprintf("+%d+%d", padding, padding),
        color = "white", size = base_font)
    }

    # J. Setup base name -------------------------------------------------------
    if (is.null(file_name)) {
      base_name <- sprintf("slice_%s_%03d", curr_plane, curr_slice)
    } else {
      base_name <- paste0(tools::file_path_sans_ext(file_name), "_", curr_plane, "_", curr_slice)
    }

    # K. Generate Alpha Mask, for later compositing ----------------------------
    if (draw_mask) {
      mask_hex <- matrix("#000000", nrow = nrow(m_slice), ncol = ncol(m_slice))
      mask_hex[m_slice != 0] <- "#FFFFFF"
      mask_img <- magick::image_read(mask_hex)
      mask_img <- magick::image_flop(magick::image_rotate(mask_img, 270))
      mask_img <- magick::image_resize(mask_img,
                                       geometry = magick::geometry_size_pixels(target_w, target_h, FALSE),
                                       filter="Catrom")
      mask_img <- magick::image_convert(mask_img, colorspace = "gray")
      magick::image_write(mask_img, path = file.path(dir_save, paste0(base_name, "_mask.png")))
    }

    # L. Save output -----------------------------------------------------------
    out_path <- file.path(dir_save, paste0(base_name, ".png"))
    output_paths <- c(output_paths, out_path)
    magick::image_write(img, path = out_path, format = "png")
  }

  # Generate the Color Bar PNG (Optional) ----------------------------------
  if (!is.null(draw_cbar)) {
    cbar_res <- 256
    pal <- grDevices::colorRampPalette(color)(cbar_res)
    label_min <- sprintf("%.2f", v_min)
    label_max <- sprintf("%.2f", v_max)
    cbar_bg <- bg_color
    target_len <- if(draw_cbar == "vertical") target_h else target_w
    strip_w  <- round(target_len * 0.5)
    strip_h  <- round(12 * scale)
    cbar_mx    <- matrix(pal, nrow = 1)
    cbar_strip <- magick::image_read(cbar_mx)
    cbar_strip <- magick::image_resize(cbar_strip,
      geometry = magick::geometry_size_pixels(width = strip_w,
                                              height = strip_h,
                                              preserve_aspect = FALSE))
    cbar_img <- magick::image_extent(cbar_strip,
                                     geometry = sprintf("%dx%d", target_len, strip_h),
                                     gravity = "Center", color = cbar_bg)
    cbar_img <- magick::image_convert(cbar_img, type = "truecoloralpha")
    cbar_img <- magick::image_annotate(cbar_img, label_min, gravity = "west",
                                       location = "+5+0", color = "white", size = base_font)
    cbar_img <- magick::image_annotate(cbar_img, label_max, gravity = "east",
                                       location = "+5+0", color = "white", size = base_font)
    if (draw_cbar == "vertical") {
      cbar_img <- magick::image_rotate(cbar_img, 270)
    }
    if (is.null(file_name)) {
      cbar_name <- "slice_cbar.png"
    } else {
      cbar_name <- sprintf("%s_cbar.png", tools::file_path_sans_ext(file_name))
    }
    cbar_out_path <- file.path(dir_save, cbar_name)
    magick::image_write(cbar_img, path = cbar_out_path, format = "png")
  }
  if (cleanup && !is.null(dir_scratch) && dir_scratch != "" && dir_scratch != getwd()) {
    unlink(dir_scratch, recursive = TRUE)
  }
  return(output_paths)
}
