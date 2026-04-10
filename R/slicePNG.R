slicePNG <- function(nii_data,
                     nii_vol = 1,
                     nii_mask = NULL,
                     mask_vol = 1,
                     plane = "coronal",
                     slice,
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
                     dir_save) {

  if (!is.null(threshold_pct) && !is.null(threshold_value)) {
    stop("Provide either threshold_pct OR threshold_value, not both.")
  }

  # 1. Setup Scratch and Local Files ------------------------------------------
  if (!dir.exists(dir_scratch)) dir.create(dir_scratch, recursive = TRUE)

  prep_file <- function(path, tag) {
    if (is.null(path) || path == "none") return(NULL)
    dest <- file.path(dir_scratch, paste0(tag, "_", gsub("\\.gz$", "", basename(path))))
    if (!file.exists(dest)) {
      if (grepl("\\.gz$", path)) {
        R.utils::gunzip(path, destname = dest, remove = FALSE, overwrite = TRUE)
      } else {
        file.copy(path, dest, overwrite = TRUE)
      }
    }
    return(dest)
  }

  local_nii <- prep_file(nii_data, "data")
  local_mask <- prep_file(nii_mask, "mask")

  # 2. Get Metadata and Thresholds --------------------------------------------
  spacing <- nifti.io::info.nii(local_nii, "spacing")
  img_vol <- nifti.io::read.nii.volume(local_nii, nii_vol)
  dims <- dim(img_vol)

  m_vol <- NULL
  if (!is.null(local_mask)) {
    m_vol <- nifti.io::read.nii.volume(local_mask, mask_vol)
    m_vol <- ((m_vol != 0) & (!is.na(img_vol)))*1
  }

  # 3. Resolve Thresholds -----------------------------------------------------
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

  # 3. Extract Slice and Calculate Physical Dimensions ------------------------
  if (plane == "coronal") {
    slice_data <- img_vol[, slice, ]
    phys_w <- dims[1] * spacing[1]; phys_h <- dims[3] * spacing[3]
  } else if (plane == "axial") {
    slice_data <- img_vol[, , slice]
    phys_w <- dims[1] * spacing[1]; phys_h <- dims[2] * spacing[2]
  } else {
    slice_data <- img_vol[slice, , ]
    phys_w <- dims[2] * spacing[2]; phys_h <- dims[3] * spacing[3]
  }

  if (is.null(local_mask)) {
    m_vol <- (slice_data[slice_data >= min(c(v_min, v_max)) & slice_data <= max(c(v_min, v_max))] != 0) * 1
  }

  # 4. Generate Main Image -----------------------------------------------------
  clamped_data <- slice_data
  clamped_data[clamped_data < v_min] <- v_min
  clamped_data[clamped_data > v_max] <- v_max
  norm_data <- (clamped_data - v_min) / (v_max - v_min)

  pal <- grDevices::colorRampPalette(color)(256)
  color_indices <- findInterval(norm_data, seq(0, 1, length.out = 256), all.inside = TRUE)

  # This is your RGB hex matrix (e.g., "#FF0000")
  hex_mx <- matrix(pal[color_indices], nrow = nrow(slice_data))
  hex_mx[is.na(hex_mx)] <- "#000000"

  #  Build the Alpha Channel from the Mask
  m_slice <- switch(plane,
                    "sagittal" = m_vol[slice, , ],
                    "coronal"  = m_vol[, slice, ],
                    "axial"    = m_vol[, , slice])
  # Create alpha values: "FF" (opaque) for data, "00" (transparent) for background
  alpha_mx <- matrix("00", nrow = nrow(slice_data), ncol = ncol(slice_data))
  alpha_mx[m_slice != 0] <- "FF"

  # Combine RGB + Alpha into a single RGBA Matrix
  # This results in hex codes like "#FF0000FF" (Opaque Red) or "#FF000000" (Transparent)
  rgba_mx <- matrix(paste0(hex_mx, alpha_mx), nrow = nrow(slice_data))

  # Create Magick Image and Composite over bg_color
  # Create the brain data with transparency ALREADY BUILT IN
  data_img <- magick::image_read(rgba_mx)
  data_img <- magick::image_flop(magick::image_rotate(data_img, 270))

  target_w <- round(phys_w * scale)
  target_h <- round(phys_h * scale)
  data_img <- magick::image_resize(data_img, geometry = magick::geometry_size_pixels(width = target_w, height = target_h, preserve_aspect = FALSE))

  # Create the final canvas and layer the pre-masked data over it
  bg_canvas <- magick::image_blank(target_w, target_h, color = bg_color)
  img <- magick::image_composite(bg_canvas, data_img, operator = "Over")


    # 5. Filenames --------------------------------------------------------------
  if (is.null(file_name)) {
    base_name <- sprintf("slice_%s_%03d", plane, slice)
  } else {
    base_name <- gsub("\\.png$", "", file_name)
  }

  # 6. Generate Alpha Mask ----------------------------------------------------
  if (draw_mask) {
    mask_hex <- matrix("#000000", nrow = nrow(m_slice), ncol = ncol(m_slice))
    mask_hex[m_slice != 0] <- "#FFFFFF"
    mask_img <- magick::image_read(mask_hex)
    mask_img <- magick::image_flop(magick::image_rotate(mask_img, 270))
    mask_img <- magick::image_resize(mask_img,
                                     geometry = magick::geometry_size_pixels(target_w, target_h, FALSE))
    # Force to grayscale to ensure CopyOpacity treats it as a single alpha channel
    mask_img <- magick::image_convert(mask_img, colorspace = "gray")
    magick::image_write(mask_img, path = file.path(dir_save, paste0(base_name, "_mask.png")))
  }

  # 7. Physical Scaling -------------------------------------------------------
  target_w <- round(phys_w * scale)
  target_h <- round(phys_h * scale)
  img <- magick::image_resize(img, geometry = magick::geometry_size_pixels(width = target_w, height = target_h, preserve_aspect = FALSE))

  # 8. Add Annotations --------------------------------------------------------
  # Calculate font sizing based on scale
  base_font <- max(10, round(10 * scale))
  padding <- 15

  if (draw_label_layer) {
    label_layer <- magick::image_blank(width = target_w, height = target_h, color = "none")
  }

  # A. Bottom Right: Scale Bar (15% width)
  if (draw_scale) {
    target_mm <- phys_w * 0.15
    scale_len_mm <- round(target_mm / 5) * 5
    if (scale_len_mm == 0) scale_len_mm <- 5
    bar_w_px <- scale_len_mm * scale
    bar_h_px <- max(1, round(0.3 * scale))
    if (draw_label_layer) {
      label_layer <- magick::image_draw(label_layer)
    } else {
      img <- magick::image_draw(img)
    }
    symbols(target_w - padding - (bar_w_px/2),
            target_h - padding - (base_font * 1.5) - (bar_h_px/2),
            rectangles = matrix(c(bar_w_px, bar_h_px), ncol=2),
            inches = FALSE, add = TRUE, bg = "white", fg = NA)
    dev.off()

    # Label UNDER the line
    if (draw_label_layer) {
      label_layer <- magick::image_annotate(label_layer, sprintf("%g mm", scale_len_mm),
                                    gravity = "southeast", color = "white", size = base_font,
                                    location = sprintf("+%d+%d", round(padding + (bar_w_px/2) - (base_font)), padding))
    } else {
      img <- magick::image_annotate(img, sprintf("%g mm", scale_len_mm),
                                    gravity = "southeast", color = "white", size = base_font,
                                    location = sprintf("+%d+%d", round(padding + (bar_w_px/2) - (base_font)), padding))
    }
  }

  # B. Bottom Left: L/R Side Label
  if (draw_side) {
    side_label <- ""
    srow_x <- unlist(nifti.io::info.nii(local_nii, "srow_x"))
    srow_y <- unlist(nifti.io::info.nii(local_nii, "srow_y"))
    srow_z <- unlist(nifti.io::info.nii(local_nii, "srow_z"))
    if (plane %in% c("coronal", "axial")) {
      side_label <- if (srow_x[1] > 0) "R" else "L"
    } else if (plane == "sagittal") {
      # For sagittal, the horizontal axis is Voxel-Y (Posterior -> Anterior)
      side_label <- if (srow_y[2] > 0) "A" else "P"
    }
    if (draw_label_layer) {
      label_layer <- magick::image_annotate(label_layer, side_label, gravity = "southwest",
                                            location = sprintf("+%d+%d", padding, padding),
                                            color = "white", size = base_font * 1.2, weight = 700)
    } else {
      img <- magick::image_annotate(img, side_label, gravity = "southwest",
                                    location = sprintf("+%d+%d", padding, padding),
                                    color = "white", size = base_font * 1.2, weight = 700)
    }
  }

  # C. Top Left: Slice Location (World mm)
  if (draw_coords) {
    # srow_x = Row 1, srow_y = Row 2, srow_z = Row 3
    hdr_field <- switch(plane, "sagittal"="srow_x", "coronal"="srow_y", "axial"="srow_z")
    tform <- unlist(nifti.io::info.nii(local_nii, hdr_field))
    # Real world mm = (index - 1) * step + offset
    world_mm <- round((slice - 1) * tform[1] + tform[4], 1)
    if (draw_label_layer) {
      label_layer <- magick::image_annotate(label_layer, sprintf("%g mm", world_mm),
                                            gravity = "northwest", location = sprintf("+%d+%d", padding, padding),
                                            color = "white", size = base_font)
    } else {
      img <- magick::image_annotate(img, sprintf("%g mm", world_mm),
                                    gravity = "northwest", location = sprintf("+%d+%d", padding, padding),
                                    color = "white", size = base_font)
    }
  }

  # 9. Save -------------------------------------------------------------------
  if (is.null(file_name)) {
    base_name <- sprintf("slice_%s_%03d", plane, slice)
  } else {
    base_name <- gsub("\\.png$", "", file_name)
  }

  out_path <- file.path(dir_save, paste0(base_name, ".png"))
  magick::image_write(img, path = out_path, format = "png")

  # 10. Generate the Color Bar PNG (Optional) ----------------------------------
  # 2. Generate the Color Bar PNG ---------------------------------------------
  if (!is.null(draw_cbar)) {
    cbar_res <- 256
    pal <- grDevices::colorRampPalette(color)(cbar_res)
    base_font <- max(10, round(10 * scale))

    label_min <- sprintf("%.2f", v_min)
    label_max <- sprintf("%.2f", v_max)
    bg_color <- color[1] # Usually black

    # Determine the "long" dimension based on the requested orientation
    # If user wants vertical, the 'length' should match the target_h
    target_len <- if(draw_cbar == "vertical") target_h else target_w

    # 1. Create a Horizontal Color Strip (50% of the target length)
    strip_w  <- round(target_len * 0.5)
    strip_h  <- round(12 * scale)
    cbar_mx    <- matrix(pal, nrow = 1)
    cbar_strip <- magick::image_read(cbar_mx)
    cbar_strip <- magick::image_resize(cbar_strip,
                                       geometry = magick::geometry_size_pixels(width = strip_w,
                                                                               height = strip_h,
                                                                               preserve_aspect = FALSE))

    # 2. Extend Canvas to full target length and add text at ends
    cbar_img <- magick::image_extent(cbar_strip,
                                     geometry = sprintf("%dx%d", target_len, strip_h),
                                     gravity = "Center", color = bg_color)

    # Labels on the horizontal version
    cbar_img <- magick::image_annotate(cbar_img, label_min, gravity = "west",
                                       location = "+5+0", color = "white", size = base_font)
    cbar_img <- magick::image_annotate(cbar_img, label_max, gravity = "east",
                                       location = "+5+0", color = "white", size = base_font)

    # 3. If vertical requested, rotate the whole thing
    if (draw_cbar == "vertical") {
      # Rotate 270 (90 deg clockwise) so Max is at the Top, Min is at the Bottom
      cbar_img <- magick::image_rotate(cbar_img, 270)
      # Note: image_rotate swaps width/height, so it now matches target_h perfectly
    }

    cbar_out_path <- file.path(dir_save, paste0(base_name, "_cbar.png"))
    magick::image_write(cbar_img, path = cbar_out_path, format = "png")
  }

  if (draw_label_layer) {
    layer_out_path <- file.path(dir_save, paste0(base_name, "_labels.png"))
    magick::image_write(label_layer, path = layer_out_path, format = "png")
  }

  return(out_path)
}
