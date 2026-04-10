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
                       plane = "coronal",
                       slice,
                       scale = 1,
                       draw_side = FALSE,
                       draw_coords = FALSE,
                       draw_scale = FALSE,
                       draw_cbar = NULL,
                       file_name = NULL,
                       dir_scratch,
                       dir_save) {

  # 1. Internal Decompression Helper ------------------------------------------
  prep_file <- function(path, tag) {
    if (is.null(path) || length(path) == 0) return(NULL)
    if (length(path) > 1) stop(sprintf("prep_file received a vector for %s. Expected single path.", tag))
    if (path == "none") return(NULL)
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

  # 2. Setup Local Files ------------------------------------------------------
  local_bg <- prep_file(bg_nii, "bg")
  if (!missing(bg_mask)) {
    local_bg_mask <- prep_file(bg_mask, "bgmask")
  } else {
    local_bg_mask <- NULL
  }

  n_fg <- 0
  if (!missing(fg_nii_list)) {
    n_fg <- length(fg_nii_list)
    local_fg <- lapply(seq_along(fg_nii_list), function(i) prep_file(fg_nii_list[i], paste0("fg", i)))
    if (!missing(fg_mask_list)) {
      local_fg_mask <- lapply(seq_along(fg_mask_list), function(i) prep_file(fg_mask_list[i], paste0("fgmask", i)))
    } else {
      local_fg_mask <- rep(NA,n_fg)
    }
  }
  if (!missing(roi_nii_list)) {
    n_rois <- length(roi_nii_list)
    local_roi_list <- lapply(seq_along(roi_nii_list), function(i) prep_file(roi_nii_list[i], paste0("roi", i)))
  }

  # set basename ---------------------------------------------------------------
  if (is.null(file_name)) {
    base_name <- sprintf("overlay_%s_%03d", plane, slice)
  } else {
    base_name <- gsub("\\.png$", "", file_name)
  }

  # 3. Resolve Background & Setup Canvas ---------------------------------------
  if (missing(bg_threshold_value)) { bg_threshold_value <- NULL }
  bg_path <- slicePNG(nii_data = local_bg, nii_mask = local_bg_mask,
                      nii_vol = bg_vol, mask_vol = bg_mask_vol,
                      plane = plane, slice = slice, scale = scale,
                      color = bg_color, bg_color = canvas_color,
                      threshold_pct = bg_threshold_pct,
                      threshold_value = bg_threshold_value,
                      draw_scale = draw_scale, draw_side = draw_side,
                      draw_coords = draw_coords, draw_label_layer = TRUE,
                      file_name = paste0(file_name, "_bg"),
                      dir_scratch = dir_scratch, dir_save = dir_scratch)
  img_stack <- magick::image_read(bg_path)
  img_stack <- magick::image_convert(img_stack, type = "truecoloralpha")

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
    for (i in 1:n_fg) {
      fg_path <- slicePNG(nii_data = local_fg[[i]],
                          nii_mask = local_fg_mask[[i]],
                          nii_vol = fg_vol_list[i],
                          mask_vol = fg_mask_vol_list[i],
                          plane = plane, slice = slice, color = fg_colors[[i]],
                          bg_color="none",
                          threshold_pct = fg_threshold_pct[[i]],
                          threshold_value = fg_threshold_value[[i]],
                          scale = scale, draw_mask = TRUE, draw_cbar = draw_cbar,
                          file_name = sprintf("%s_fg%d", file_name, i),
                          dir_scratch = dir_scratch, dir_save = dir_save)
      fg_img <- magick::image_read(fg_path)
      fg_img <- magick::image_convert(fg_img, type = "truecoloralpha")
      fg_mask <- magick::image_read(gsub("\\.png$", "_mask.png", fg_path))
      fg_mask <- magick::image_convert(fg_mask, colorspace = "gray")
      fg_trans <- magick::image_composite(fg_img, fg_mask, operator = "CopyOpacity")
      if (fg_alphas[i] < 1) {
        fg_trans <- magick::image_colorize(fg_trans, opacity = (1 - fg_alphas[i]) * 100, color = "transparent")
      }
      img_stack <- magick::image_composite(img_stack, fg_trans, operator = "Over")
      file.remove(fg_path, gsub("\\.png$", "_mask.png", fg_path))
    }
  }

  # 5. Multiple ROI Outlines --------------------------------------------------
  if (!missing(roi_nii_list)) {
    if (length(roi_vol_list) == 1) { roi_vol_list <- rep(roi_vol_list, n_rois) }
    roi_colors <- if (length(roi_color) == 1) rep(roi_color, n_rois) else roi_color
    for (r_idx in 1:n_rois) {
      roi_vol <- nifti.io::read.nii.volume(prep_file(roi_nii_list[[r_idx]], "roi"), roi_vol_list[r_idx])
      if (roi_value != "all") roi_vol <- (roi_vol == as.numeric(roi_value)) * 1
      r_slice <- switch(plane, "coronal"=roi_vol[,slice,], "axial"=roi_vol[,,slice], "sagittal"=roi_vol[slice,,])
      interior <- (r_slice > 0) & (cbind(r_slice[,-1], 0) > 0) & (cbind(0, r_slice[,-ncol(r_slice)]) > 0) &
        (rbind(r_slice[-1,], 0) > 0) & (rbind(0, r_slice[-nrow(r_slice),]) > 0)
      outline_mx <- (r_slice > 0 & !interior)
      if (sum(outline_mx) > 0) {
        roi_hex <- matrix("transparent", nrow = nrow(r_slice), ncol = ncol(r_slice))
        roi_hex[outline_mx] <- roi_colors[r_idx]
        roi_img <- magick::image_flop(magick::image_rotate(magick::image_read(roi_hex), 270))
        info <- magick::image_info(img_stack)
        roi_img <- magick::image_resize(roi_img, geometry = magick::geometry_size_pixels(width = info$width, height = info$height, preserve_aspect = FALSE))
        img_stack <- magick::image_composite(img_stack, roi_img, operator = "Over")
      }
    }
  }

  # 6. Final Annotations ------------------------------------------------------
  label_layer_path <- file.path(dir_scratch, paste0(base_name, "_bg_labels.png"))
  if (file.exists(label_layer_path)) {
    label_img <- magick::image_read(label_layer_path)
    # Composite the labels (White text/scale bar) on top of the entire stack
    img_stack <- magick::image_composite(img_stack, label_img, operator = "Over")
    # Cleanup the label layer from scratch
    file.remove(label_layer_path)
  }

  # 7. Output -----------------------------------------------------------------
  if (is.null(file_name)) {
    base_name <- sprintf("overlay_%s_%03d", plane, slice)
  } else {
    base_name <- gsub("\\.png$", "", file_name)
  }
  out_path <- file.path(dir_save, paste0(base_name, ".png"))
  magick::image_write(img_stack, path = out_path, format = "png")

  if (!is.null(draw_cbar)) {
    cbar_files <- list.files(dir_save, pattern = "_cbar.png$", full.names = TRUE)
  }

  return(out_path)
}
