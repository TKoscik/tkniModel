drawEffect <- function(model,
                       effect_name,
                       axes = NULL,
                       labels = NULL,
                       plot_colors = get.color(),
                       save_dir,
                       file_name,
                       img_format = "png",
                       img_w = 10,
                       img_unit = "cm",
                       img_dpi = 600,
                       save_plot = TRUE,
                       return_plot = FALSE) {

  # Ensure the save directory exists
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # 1. Model Parsing and Random Effect Adjustment ----------------------------
  ## Note on Partialling Out Random Effects: For mixed-effects models (lmer),
  ## this function visualizes the relationship after adjusting for the
  ## hierarchical structure of the data. It calculates the "partialled-out"
  ## dependent variable by subtracting the Best Linear Unbiased Predictors
  ## (BLUPs) of the random effects from the observed values. This removes the
  ## variance attributed to specific subjects or groups, effectively showing
  ## the "fixed effect" relationship—or what the data would look like if every
  ## group had the same intercept and slope. This ensures the visual trend
  ## matches the statistical coefficients reported by the model.
  if (inherits(model, "merModLmerTest") || inherits(model, "lmerMod")) {
    tempf <- model@frame
    # Extract Dependent Variable (DV) name from the formula
    dv <- as.character(formula(model))[2]

    tempf[[paste0(dv, ".dummy")]] <- tempf[[dv]]
    temp_ranef <- lme4::ranef(model)

    # Iterate through random effects to subtract them from the DV
    for (i in 1:length(temp_ranef)) {
      which_vars <- names(temp_ranef)[[i]]
      which_vars <- unlist(strsplit(which_vars, split = ":"))

      # Create a dummy ID column for the random effect grouping
      tempf[[paste0("dummy", i)]] <- do.call(paste, c(tempf[which_vars], sep = ":"))

      # Match the random effect value to the row and subtract
      # This effectively 'corrects' the DV for that random level
      re_values <- temp_ranef[[i]][[1]] # Assuming intercept-based ranef
      names(re_values) <- rownames(temp_ranef[[i]])

      tempf[[paste0(dv, ".dummy")]] <- tempf[[paste0(dv, ".dummy")]] -
        re_values[as.character(tempf[[paste0("dummy", i)]])]
    }

    # Update the DV in our local frame and re-run the model locally for plotting
    tempf[[dv]] <- tempf[[paste0(dv, ".dummy")]]
    model <- stats::update(model, data = tempf)

  } else {
    # Standard lm/glm
    tempf <- model$model
  }

  # 2. Determine Plotting Mapping --------------------------------------------
  if (is.null(axes)) {
    # Automatic mapping based on the interaction order
    ef_vars <- unlist(strsplit(effect_name, split = ":"))
    ef_num <- length(ef_vars)

    # Check variable types (numeric vs factor)
    var_type <- sapply(tempf[ef_vars], function(x) class(x)[1])

    # Map variables to axes: c(x, y, color, facet1, facet2)
    axes <- switch(as.character(ef_num),
                   "1" = c(ef_vars[1], "fit"),
                   "2" = c(ef_vars[1], "fit", ef_vars[2]),
                   "3" = c(ef_vars[1], "fit", ef_vars[2], ef_vars[3]),
                   "4" = c(ef_vars[1], "fit", ef_vars[2], ef_vars[3], ef_vars[4]),
                   stop("Cannot plot interactions larger than 4-way"))
  } else {
    # User provided custom axes
    ef_vars <- axes[-which(axes == "fit")]
    ef_num <- length(ef_vars)
    var_type <- sapply(tempf[ef_vars], function(x) class(x)[1])
  }

  # 3. Set Labels -------------------------------------------------------------
  if (is.null(labels)) {
    # Capitalize the first letter of each axis name
    labels <- tools::toTitleCase(axes)
  }

  # 4. Determine Effect Levels for Prediction ---------------------------------
  effect_levels <- numeric(ef_num)
  for (i in 1:ef_num) {
    current_var <- ef_vars[i]
    current_class <- var_type[[i]]

    if (i == 1) {
      # Primary X-axis
      if (current_class == "factor") {
        effect_levels[i] <- nlevels(tempf[[current_var]])
      } else if (current_class %in% c("numeric", "integer")) {
        effect_levels[i] <- 100 # High resolution for smooth lines
      } else {
        stop(sprintf("Cannot parse class of variable: %s", current_var))
      }
    } else {
      # Interaction/Grouping variables
      if (current_class == "factor") {
        effect_levels[i] <- nlevels(tempf[[current_var]])
      } else if (current_class %in% c("numeric", "integer")) {
        # If numeric, pick 3 points (typically -1SD, Mean, +1SD logic) or unique values
        n_unique <- length(unique(tempf[[current_var]]))
        effect_levels[i] <- ifelse(n_unique >= 3, 3, n_unique)
      } else {
        stop(sprintf("Cannot parse class of variable: %s", current_var))
      }
    }
  }

  # 5. Calculate Predicted Effects --------------------------------------------
  # We use the 'effects' package to generate the predicted means and CIs
  # We construct the xlevels list dynamically
  xlevs <- list()
  for (i in 1:ef_num) {
    xlevs[[ef_vars[i]]] <- effect_levels[i]
  }

  # Generate the effect data frame
  # This replaces the switch(as.character(ef.num)...) and eval(parse()) logic
  ef_obj <- effects::effect(effect_name, model, xlevels = xlevs)
  ef <- as.data.frame(ef_obj)

  # 6. Formatting and Labeling for Plotting -----------------------------------
  # Apply rounding to numeric grouping variables for cleaner legend labels
  if (length(axes) > 2) {
    group_var <- axes[3]
    if (is.numeric(ef[[group_var]])) {
      ef[[group_var]] <- signif(ef[[group_var]], digits = 3)
    }
  }

  # Update column names for facet variables to match user-provided labels
  # This ensures the plot headers (e.g., in facet_wrap) look professional
  if (ef_num >= 3) {
    # axes[4] is the 1st facet variable
    col_idx4 <- which(colnames(ef) == axes[4])
    if(length(col_idx4) > 0) colnames(ef)[col_idx4] <- labels[4]

    if (ef_num == 4) {
      # axes[5] is the 2nd facet variable
      col_idx5 <- which(colnames(ef) == axes[5])
      if(length(col_idx5) > 0) colnames(ef)[col_idx5] <- labels[5]
    }
  }

  # 7. Plotting Categorical Effects (X is Factor) -----------------------------
  if (var_type[1] == "factor") {

    # Base plot common to all interaction levels
    p <- ggplot2::ggplot(ef, ggplot2::aes_string(x = axes[1], y = axes[2],
                                                 ymin = "lower", ymax = "upper")) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = labels[1], y = labels[2]) +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                     axis.text  = ggplot2::element_text(size = 12),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
                     legend.position = "bottom",
                     legend.title = ggplot2::element_text(size = 14),
                     legend.text = ggplot2::element_text(size = 12))

    # Add layers based on complexity (1-way to 4-way)
    if (ef_num == 1) {
      p <- p + ggplot2::geom_pointrange(size = 1, shape = 18)

    } else {
      # For 2, 3, and 4-way, we add color/fill and dodging
      p <- p + ggplot2::aes_string(fill = paste0("factor(", axes[3], ")"),
                                   color = paste0("factor(", axes[3], ")")) +
        ggplot2::geom_pointrange(size = 1, shape = 18,
                                 position = ggplot2::position_dodge(width = 0.5)) +
        ggplot2::scale_fill_manual(values = plot_colors, name = labels[3]) +
        ggplot2::scale_color_manual(values = plot_colors, name = labels[3])

      if (ef_num == 3) {
        # 3-way: Facet Grid by the 4th variable (labeled in step 6)
        p <- p + ggplot2::facet_grid(stats::as.formula(paste0(". ~ `", labels[4], "`")),
                                     labeller = ggplot2::label_both) +
          ggplot2::theme(strip.background = ggplot2::element_rect(fill = '#ffffff', size = 0.5),
                         strip.text = ggplot2::element_text(size = 12))

      } else if (ef_num == 4) {
        # 4-way: Grid by 5th (rows) and 4th (cols) variables
        p <- p + ggplot2::facet_grid(stats::as.formula(paste0("`", labels[5], "` ~ `", labels[4], "`")),
                                     labeller = ggplot2::label_both) +
          ggplot2::theme(strip.background = ggplot2::element_rect(fill = '#ffffff', size = 0.5),
                         strip.text = ggplot2::element_text(size = 12))
      }
    }
  }  # 8. Plotting Numeric Effects (X is Continuous) ----------------------------
  else if (var_type[1] %in% c("numeric", "integer")) {

    # Base plot for continuous X
    p <- ggplot2::ggplot(ef, ggplot2::aes_string(x = axes[1], y = axes[2],
                                                 ymin = "lower", ymax = "upper")) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = labels[1], y = labels[2]) +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                     axis.text  = ggplot2::element_text(size = 12),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
                     legend.position = "bottom",
                     legend.title = ggplot2::element_text(size = 14),
                     legend.text = ggplot2::element_text(size = 12))

    if (ef_num == 1) {
      p <- p +
        ggplot2::geom_ribbon(alpha = 0.25, color = "transparent") +
        ggplot2::geom_line(linewidth = 1)

    } else {
      # For 2, 3, and 4-way, apply color/fill to the ribbon and line
      p <- p + ggplot2::aes_string(fill = paste0("factor(", axes[3], ")"),
                                   color = paste0("factor(", axes[3], ")")) +
        ggplot2::geom_ribbon(color = "transparent", alpha = 0.25) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_fill_manual(values = plot_colors, name = labels[3]) +
        ggplot2::scale_color_manual(values = plot_colors, name = labels[3])

      if (ef_num == 3) {
        p <- p + ggplot2::facet_grid(stats::as.formula(paste0(". ~ `", labels[4], "`")),
                                     labeller = ggplot2::label_both) +
          ggplot2::theme(strip.background = ggplot2::element_rect(fill = '#ffffff', linewidth = 0.5),
                         strip.text = ggplot2::element_text(size = 12))

      } else if (ef_num == 4) {
        p <- p + ggplot2::facet_grid(stats::as.formula(paste0("`", labels[5], "` ~ `", labels[4], "`")),
                                     labeller = ggplot2::label_both) +
          ggplot2::theme(strip.background = ggplot2::element_rect(fill = '#ffffff', linewidth = 0.5),
                         strip.text = ggplot2::element_text(size = 12))
      }
    }
  }

  # 9. Dynamic Height Calculation and Saving ---------------------------------
  if (save_plot) {
    # Calculate height based on interaction complexity and number of facets
    img_h <- switch(as.character(ef_num),
                    "1" = (img_w / 6) * 5,
                    "2" = (img_w / 6) * 6,
                    "3" = (img_w / (3 * effect_levels[3])) * 6,
                    "4" = (img_w / (3 * effect_levels[4])) * (6 * effect_levels[3])
    )

    out_file <- file.path(save_dir, paste0(file_name, ".", img_format))

    ggplot2::ggsave(
      filename = out_file,
      plot     = p,
      width    = img_w,
      height   = img_h,
      units    = img_unit,
      dpi      = img_dpi
    )

    if (verbose) message(sprintf("Effect plot saved to: %s", out_file))
  }

  # 10. Returns ---------------------------------------------------------------
  if (return_plot) {
    return(p)
  } else {
    return(invisible(file.path(save_dir, paste0(file_name, ".", img_format))))
  }
}
