drawOverlay <- function(anat_nii,
                        over_nii, over_vol, over_color,
                        mask_nii = "none",
                        mask_vol = 1,
                        roi_nii = NULL,
                        roi_val = NULL,
                        roi_color = "#ff64ff",
                        orientation = "coronal",
                        save_dir = getwd(),
                        file_name = "overlay",
                        img_format = "png",
                        img_w = NULL,
                        img_unit = "cm",
                        img_dpi = 600,
                        save_plot = TRUE,
                        return_plot = FALSE) {

  # 1. Normalize Input Lengths ------------------------------------------------
  n_check <- c(length(over_nii), length(over_color), length(mask_nii))
  n_overlays <- max(n_check)

  if (!all(n_check == n_overlays)) {
    # Expand over_nii and over_vol
    if (length(over_nii) != n_overlays) {
      temp <- over_nii
      over_nii <- vector("list", n_overlays)
      for (i in 1:n_overlays) { over_nii[[i]] <- temp }
      over_vol <- rep(over_vol, n_overlays)
    }

    # Expand over_color
    if (length(over_color) != n_overlays) {
      temp <- over_color
      over_color <- vector("list", n_overlays)
      for (i in 1:n_overlays) { over_color[[i]] <- temp }
    }

    # Expand mask_nii and mask_vol
    if (length(mask_nii) != n_overlays) {
      temp_mask <- mask_nii
      mask_nii <- vector("list", n_overlays)
      for (i in 1:n_overlays) { mask_nii[[i]] <- temp_mask }

      temp_vol <- mask_vol
      mask_vol <- vector("list", n_overlays)
      for (i in 1:n_overlays) { mask_vol[[i]] <- temp_vol }
    }
  }

  if (missing(img_format)) { img_format <- "png" }

  # Load Anatomical ----
  img.anat <- read.nii.volume(anat.nii, 1)
  img.anat[img.anat==0] <- NA
  img.dims <- dim(img.anat)

  # Load Overlays ----
  img.over <- vector("list", length=n.overlays)
  for (i in 1:n.overlays) { img.over[[i]] <- read.nii.volume(over.nii[[i]], over.vol[i]) }

  # Load mask ----
  img.mask <- vector("list", length=n.overlays)
  for (i in 1:n.overlays) {
    if (mask.nii[[i]] == "none") {
      img.mask[[i]] <- array(as.numeric(img.over[[i]] != 0), dim=dim(img.anat))
    } else {
      if (mask.vol[[i]] == "all") {
        mask.vols <- 1:(nii.dims(mask.nii[[i]])[4])
      } else if (is.numeric(mask.vol[[i]])) {
        mask.vols <- mask.vol[[i]]
      } else { stop("Cannot parse mask volumes") }

      img.mask[[i]] <- array(0, dim=dim(img.anat))
      for (j in mask.vols) {
        img.mask[[i]] <- img.mask[[i]] + read.nii.volume(mask.nii[[i]], j)
      }
      img.mask[[i]][img.mask[[i]]>1] <- 1
    }
    img.over[[i]][img.mask[[i]]==0] <- NA
  }

  # Load ROIs ----
  if (!is.null(roi.nii)) {
    img.roi <- read.nii.volume(roi.nii, 1)
    ex.rois <- c(0, which(!(1:max(img.roi) %in% roi.val)))
    for (i in ex.rois) { img.roi[img.roi==i] <- NA }
  }

  # 3. Slice Selection Logic --------------------------------------------------
  slices <- numeric(0)
  if (!is.null(roi_nii)) {
    all_mask <- !is.na(img_roi)
  } else {
    all_mask <- array(0, dim = dim(img_over[[1]]))
  }

  for (i in 1:n_overlays) {
    all_mask <- all_mask + !is.na(img_over[[i]])
  }

  # Normalize orientation string
  orientation <- tolower(orientation)

  # Determine which plane to scan
  scan_dim <- switch(orientation, "coronal" = 2, "c" = 2, "axial" = 3, "a" = 3, "sagittal" = 1, "s" = 1, stop("Invalid orientation"))

  # Find slices with data
  for (i in 1:dim(all_mask)[scan_dim]) {
    slice_data <- if(scan_dim == 1) all_mask[i, , ] else if(scan_dim == 2) all_mask[, i, ] else all_mask[, , i]
    if (sum(slice_data, na.rm = TRUE) > 0) {
      slices <- c(slices, i)
    }
  }

  # Trim edges and calculate grid size
  if (length(slices) > 2) slices <- slices[-c(1, length(slices))]

  possible_lengths <- c(Inf, 40, 35, 30, 24, 20, 15, 12, 6, 4, 3, 2, 1)
  slice_length <- possible_lengths[min(which(length(slices) >= possible_lengths))]

  n_row <- switch(as.character(slice_length), "40"=5, "35"=5, "30"=5, "24"=4, "20"=4, "15"=3, "12"=3, "6"=2, "4"=2, "3"=3, "2"=2, "1"=1)
  n_col <- switch(as.character(slice_length), "40"=8, "35"=7, "30"=6, "24"=6, "20"=5, "15"=5, "12"=4, "6"=3, "4"=2, "3"=1, "2"=1, "1"=1)

  # Evenly space the slices to fit the grid
  slices <- slices[round(seq(from = 1, to = length(slices), length.out = slice_length))]

  # 4. Generate Raster Montage ------------------------------------------------
  # Helper for the internal melt/expand.grid
  simple_melt <- function(mx) {
    df <- expand.grid(Var1 = 1:nrow(mx), Var2 = 1:ncol(mx))
    df$value <- as.vector(mx)
    return(df)
  }

  if (orientation %in% c("coronal", "c")) {
    img_idx <- simple_melt(img_anat[floor(dim(img_anat)[1]/2), , ])
    img_anat_raster <- slicesToRaster(img_anat[, slices, ], n_row, n_col, "coronal")
    for (i in 1:n_overlays) img_over[[i]] <- slicesToRaster(img_over[[i]][, slices, ], n_row, n_col, "coronal")
    if (!is.null(roi_nii)) img_roi_raster <- slicesToRaster(img_roi[, slices, ], n_row, n_col, "coronal")

  } else if (orientation %in% c("axial", "a")) {
    img_idx <- simple_melt(img_anat[, floor(dim(img_anat)[2]/2), ])
    img_anat_raster <- slicesToRaster(img_anat[, , slices], n_row, n_col, "axial")
    for (i in 1:n_overlays) img_over[[i]] <- slicesToRaster(img_over[[i]][, , slices], n_row, n_col, "axial")
    if (!is.null(roi_nii)) img_roi_raster <- slicesToRaster(img_roi[, , slices], n_row, n_col, "axial")

  } else if (orientation %in% c("sagittal", "s")) {
    img_idx <- simple_melt(img_anat[, , floor(dim(img_anat)[3]/2)])
    img_idx <- simple_melt(img_idx) # Match double melt in original if needed
    img_anat_raster <- slicesToRaster(img_anat[slices, , ], n_row, n_col, "sagittal")
    for (i in 1:n_overlays) img_over[[i]] <- slicesToRaster(img_over[[i]][slices, , ], n_row, n_col, "sagittal")
    if (!is.null(roi_nii)) img_roi_raster <- slicesToRaster(img_roi[slices, , ], n_row, n_col, "sagittal")
  }

  # 5. Calculate Real-World Coordinates --------------------------------------
  # Setup the label dataframe structure for the montage
  world_labels <- data.frame(
    yvar = sort(rep(seq(0, 1 - 1/n_row, length.out = n_row), n_col), decreasing = TRUE),
    xvar = rep(seq(0, 1 - 1/n_col, length.out = n_col), n_row),
    labels = numeric(n_row * n_col)
  )

  # Get dimensions for scaling the label positions on the raster plot
  anat_dims <- nifti.io::nii.dims(local_anat)

  if (orientation %in% c("coronal", "c")) {
    tform <- unlist(nifti.io::nii.hdr(local_anat, "srow_y"))
    world_labels$xvar <- world_labels$xvar * anat_dims[1] * n_col + anat_dims[1] / 2
    world_labels$yvar <- world_labels$yvar * anat_dims[3] * n_row - 3
    # Coordinate = (0-based Index * Scale) + Offset
    world_labels$labels <- (slices - 1) * tform[2] + tform[4]

  } else if (orientation %in% c("axial", "a")) {
    tform <- unlist(nifti.io::nii.hdr(local_anat, "srow_z"))
    world_labels$xvar <- world_labels$xvar * anat_dims[1] * n_col + anat_dims[1] / 2
    world_labels$yvar <- world_labels$yvar * anat_dims[2] * n_row - 3
    world_labels$labels <- (slices - 1) * tform[3] + tform[4]

  } else if (orientation %in% c("sagittal", "s")) {
    tform <- unlist(nifti.io::nii.hdr(local_anat, "srow_x"))
    world_labels$xvar <- world_labels$xvar * anat_dims[2] * n_col + anat_dims[2] / 2
    world_labels$yvar <- world_labels$yvar * anat_dims[3] * n_row - 3
    world_labels$labels <- (slices - 1) * tform[1] + tform[4]
  }

  # Round labels for clean display in the figure
  world_labels$labels <- round(world_labels$labels, 0)

  # 6. Map Overlay Values to RGB Colors --------------------------------------
  for (i in 1:n_overlays) {
    # Generate the color interpolation function
    col_func <- grDevices::colorRamp(over_color[[i]])

    # Normalize values between 0 and 1 for the color ramp
    v_min <- min(img_over[[i]]$value, na.rm = TRUE)
    v_max <- max(img_over[[i]]$value, na.rm = TRUE)

    # Avoid division by zero if max == min
    if (v_max == v_min) {
      img_over[[i]]$value_scale <- 0
    } else {
      img_over[[i]]$value_scale <- (img_over[[i]]$value - v_min) / (v_max - v_min)
    }

    # Find non-NA indices
    valid_idx <- which(!is.na(img_over[[i]]$value_scale))

    # Initialize color columns
    img_over[[i]]$red   <- 0
    img_over[[i]]$green <- 0
    img_over[[i]]$blue  <- 0

    if (length(valid_idx) > 0) {
      # Map the scaled values to RGB [0-255] then scale to [0-1]
      rgb_vals <- col_func(img_over[[i]]$value_scale[valid_idx]) / 255

      img_over[[i]]$red[valid_idx]   <- rgb_vals[, 1]
      img_over[[i]]$green[valid_idx] <- rgb_vals[, 2]
      img_over[[i]]$blue[valid_idx]  <- rgb_vals[, 3]
    }

    # Clamp values to [0, 1] and handle NAs
    img_over[[i]]$red   <- pmin(pmax(img_over[[i]]$red, 0, na.rm = TRUE), 1)
    img_over[[i]]$green <- pmin(pmax(img_over[[i]]$green, 0, na.rm = TRUE), 1)
    img_over[[i]]$blue  <- pmin(pmax(img_over[[i]]$blue, 0, na.rm = TRUE), 1)

    # Calculate alpha based on whether there is any color present
    img_over[[i]]$alpha <- as.numeric((img_over[[i]]$red + img_over[[i]]$green + img_over[[i]]$blue) > 0)
  }

  # 7. Merge RGB Layers -------------------------------------------------------
  plotf <- img_anat_raster
  plotf$r <- 0; plotf$g <- 0; plotf$b <- 0; plotf$a <- 0

  for (i in 1:n_overlays) {
    plotf$r <- plotf$r + img_over[[i]]$red
    plotf$g <- plotf$g + img_over[[i]]$green
    plotf$b <- plotf$b + img_over[[i]]$blue
    plotf$a <- plotf$a + img_over[[i]]$alpha
  }

  # Clamp alpha to 1 to prevent "super-opaque" voxels
  plotf$a[plotf$a > 1] <- 1

  # 8. Handle ROI Polygons ----------------------------------------------------
  if (!is.null(roi_nii)) {
    # Using 'raster' and 'sp' logic as in your original snippet
    # This creates the clean outlines around the ROI
    roi_raster <- raster::rasterFromXYZ(img_roi_raster)
    roi_poly <- ggplot2::fortify(raster::rasterToPolygons(roi_raster, dissolve = TRUE))
  }

  # 9. Main Plot Construction -------------------------------------------------
  plot_img <- ggplot2::ggplot(plotf, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    # Layer 1: Anatomical grayscale background
    ggplot2::geom_raster() +
    # Layer 2: Colored statistical overlays
    ggplot2::geom_raster(inherit.aes = FALSE,
                         data = plotf,
                         ggplot2::aes(x = Var1, y = Var2),
                         fill = grDevices::rgb(red = plotf$r,
                                               green = plotf$g,
                                               blue = plotf$b,
                                               alpha = plotf$a),
                         na.rm = TRUE) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient(low = "#000000", high = "#ffffff", na.value = "transparent") +
    # Theme settings to remove all axes, ticks, and backgrounds
    ggplot2::theme(legend.position = "none",
                   legend.spacing = ggplot2::unit(0, "null"),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   plot.margin = ggplot2::unit(c(0, 0, 0, 0), "null"),
                   plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   panel.spacing = ggplot2::unit(c(0, 0, 0, 0), "null"),
                   panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank()) +
    # Add coordinate labels
    ggplot2::geom_text(inherit.aes = FALSE,
                       data = world_labels,
                       ggplot2::aes(x = xvar, y = yvar, label = labels),
                       color = "#484848", size = 3)

  # 10. Add ROI Outlines -----------------------------------------------------
  if (!is.null(roi_nii)) {
    plot_img <- plot_img +
      ggplot2::geom_path(inherit.aes = FALSE, data = roi_poly,
                         ggplot2::aes(x = long, y = lat, group = group),
                         linewidth = 0.25, alpha = 0.5,
                         color = roi_color, linetype = "solid")
  }

  # 11. Create Index/Locator Plot --------------------------------------------
  plot_idx <- ggplot2::ggplot(img_idx, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient(low = "#000000", high = "#ffffff", na.value = "transparent") +
    ggplot2::theme(legend.position = "none",
                   legend.spacing = ggplot2::unit(0, "null"),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   plot.margin = ggplot2::unit(c(0, 0, 0, 0), "null"),
                   plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   panel.spacing = ggplot2::unit(c(0, 0, 0, 0), "null"),
                   panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 14)) +
    ggplot2::annotate("text", x = -1, y = -2, label = 1, size = n_row * 0.48)

  # Orientation-specific blue slice lines and dimension markers
  if (orientation %in% c("coronal", "c")) {
    plot_idx <- plot_idx +
      ggplot2::annotate("segment", x = slices, xend = slices, y = 0, yend = Inf, color = "#0000ff", linewidth = 0.25) +
      ggplot2::annotate("text", x = img_dims[2] - 4, y = -2, label = img_dims[2], size = n_row * 0.48) +
      ggplot2::annotate("segment", x = 5, xend = img_dims[2] - 15, y = -2, yend = -2,
                        arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "npc")))

  } else if (orientation %in% c("axial", "a")) {
    plot_idx <- plot_idx +
      ggplot2::annotate("segment", x = 0, xend = Inf, y = slices, yend = slices, color = "#0000ff", linewidth = 0.25) +
      ggplot2::annotate("text", y = img_dims[1] - 4, x = -2, label = img_dims[1], size = n_row * 0.48) +
      ggplot2::annotate("segment", y = 5, yend = img_dims[1] - 15, x = -2, xend = -2,
                        arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "npc")))

  } else if (orientation %in% c("sagittal", "s")) {
    plot_idx <- plot_idx +
      ggplot2::annotate("segment", x = slices, xend = slices, y = 0, yend = Inf, color = "#0000ff", linewidth = 0.25) +
      ggplot2::annotate("text", x = img_dims[3] - 4, y = -2, label = img_dims[3], size = n_row * 0.48) +
      ggplot2::annotate("segment", x = 5, xend = img_dims[3] - 15, y = -2, yend = -2,
                        arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "npc")))
  }

  # 12. Create Color Bars ----------------------------------------------------
  plot_cbar <- vector("list", n_overlays)
  scale_res <- 500
  img_cbar <- data.frame(x = 1, y = 1:scale_res)

  for (i in 1:n_overlays) {
    # Extract min/max from the rasterized overlay values
    v_max <- max(img_over[[i]]$value, na.rm = TRUE)
    v_min <- min(img_over[[i]]$value, na.rm = TRUE)

    cbar_labels <- data.frame(
      xval = c(1, 1),
      yval = c(scale_res + 0.2 * scale_res, -0.2 * scale_res),
      the_labels = c(round(v_max, 3), round(v_min, 3))
    )

    plot_cbar[[i]] <- ggplot2::ggplot(img_cbar, ggplot2::aes(x = x, y = y, fill = y)) +
      ggplot2::geom_raster() +
      # Ratio maintains the vertical "strip" look regardless of the montage size
      ggplot2::coord_equal(ratio = ((n_col + 1) / 8) * 0.05, expand = FALSE) +
      ggplot2::scale_fill_gradientn(colors = over_color[[i]]) +
      ggplot2::geom_text(inherit.aes = FALSE, data = cbar_labels,
                         ggplot2::aes(x = xval, y = yval, label = the_labels),
                         angle = 90, size = n_row * 0.5) +
      ggplot2::ylim(c(-0.4 * scale_res, scale_res + 0.4 * scale_res)) +
      ggplot2::theme(legend.position = "none",
                     legend.spacing = ggplot2::unit(0, "null"),
                     axis.title = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     plot.margin = ggplot2::unit(c(0, 0, 0, 0), "null"),
                     plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                     panel.spacing = ggplot2::unit(c(0, 0, 0, 0), "null"),
                     panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                     panel.grid = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 14))
  }

  # 13. Assemble the Final Layout ---------------------------------------------
  # Construct the layout matrix:
  # Column 1-N: Main montage (ID 1)
  # Last Column: Stacked colorbars (IDs 2 to N+1) and the Index plot (ID N+2)
  lay_main <- matrix(1, nrow = n_row * n_overlays, ncol = n_col * n_overlays)

  # Stack colorbars vertically in the side panel
  cbar_rows <- n_overlays * (n_row - 1)
  lay_side_cbars <- matrix(sort(rep(2:(n_overlays + 1), cbar_rows / n_overlays)),
                           nrow = cbar_rows, ncol = n_overlays)

  # Place Index plot at the bottom of the side panel
  lay_side_idx <- matrix(n_overlays + 2, nrow = n_overlays, ncol = n_overlays)

  lay <- cbind(lay_main, rbind(lay_side_cbars, lay_side_idx))

  # Combine all plot objects into a single list for arrangeGrob
  all_plots <- list(plot_img)
  for (i in 1:n_overlays) { all_plots[[i+1]] <- plot_cbar[[i]] }
  all_plots[[n_overlays + 2]] <- plot_idx

  # Assemble the grob
  plot_all <- do.call(gridExtra::arrangeGrob, c(all_plots, list(layout_matrix = lay)))

  # 14. Determine Image Dimensions --------------------------------------------
  if (is.null(img_w)) {
    img_w <- (n_col + 1) * 2
    img_h <- (n_row) * 2
  } else {
    img_h <- (n_row * img_w) / (n_col + 1)
  }

  # 15. Save and Return -------------------------------------------------------
  if (save_plot) {
    out_file <- file.path(save_dir, paste0(file_name, ".", img_format))
    ggplot2::ggsave(filename = out_file, plot = plot_all,
                    width = img_w, height = img_h, units = img_unit, dpi = img_dpi)
    if (verbose) message(sprintf("Overlay saved to: %s", out_file))
  }

  # Clean up temp unzipped files before exiting
  unlink(plot_scratch, recursive = TRUE)

  if (return_plot) {
    return(plot_all)
  } else {
    return(invisible(file.path(save_dir, paste0(file_name, ".", img_format))))
  }
}
