modelVoxel <- function(nii_data,
                       df_data,
                       id_var="participant_id",
                       var_factor,
                       model_pfx,
                       model_fcn,
                       model_libraries = NULL,
                       roi_nii,
                       dir_save,
                       dir_scratch=NULL,
                       restart_log=TRUE,
                       rand_order=TRUE,
                       num_cores=parallel::detectCores()-1,
                       verbose=TRUE,
                       debug=NA,
                       cleanup=TRUE) {

  # check for missing data -----------------------------------------------------
  if (verbose) { message("Checking required inputs") }
  missing_input <- FALSE
  if (missing(nii_data)) { message("Error: 'nii_data' is required."); missing_input <- TRUE }
  if (missing(df_data)) { message("Error: 'df_data' is required."); missing_input <- TRUE }
  if (missing(model_fcn)) { message("Error: 'model_fcn' is required."); missing_input <- TRUE }
  if (missing_input) {
    stop("Missing required arguments. Please check the messages above and try again.", call. = FALSE)
  }

  if (missing(model_pfx)) { model_pfx <- sprintf("model_%s", format(Sys.time(), "%Y%m%dT%H%M%S"))}

  # setup scratch directory ------------------------------------------------------
  if (verbose) { message("Setting scratch directory") }
  if (is.null(dir_scratch)) { dir_scratch <- file.path(dir_save, "scratch") }
  if (!dir.exists(dir_scratch)) { dir.create(dir_scratch, recursive = TRUE) }
  if (verbose) { message(sprintf("Work directory set to: %s", dir_scratch)) }

  # load required libraries ------------------------------------------------------
  if (verbose) { message("Loading libraries") }
  core_libs <- c("doParallel", "nifti.io", "tools", "R.utils")
  all_libs <- unique(c(core_libs, model_libraries))
  for (lib in all_libs) {
    if (!require(lib, character.only = TRUE, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed. Please install it to continue.", lib))
    }
  }

  # load data frame for analysis -------------------------------------------------
  if (verbose) { message(sprintf("Reading data frame: %s", df_data)) }
  ext <- file_ext(df_data)
  if (ext == "csv") {
    pf <- read.csv(df_data, stringsAsFactors = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    pf <- read.delim(df_data, stringsAsFactors = FALSE)
  } else {
    pf <- read.table(df_data, header = TRUE, stringsAsFactors = FALSE)
  }

  ## make sure IDs are factors, and set the order of groups if not alphabetical
  if (verbose) { message("Converting IDs to factors") }
  for (id in id_var) {
    if (id %in% names(pf)) {
      pf[[id]] <- as.factor(pf[[id]])
    } else {
      warning(sprintf("id_var '%s' not found in data frame.", id))
    }
  }

  ## Setup ordered factors (ordered/specific levels)
  # Expected format: c("variable:levels,labels,", etc...)
  if (!missing(var_factor) && !is.null(var_factor)) {
    if (verbose) { message("Setting up variables as factors") }
    for (item in var_factor) {
      parts <- unlist(strsplit(item, ":"))
      var_name <- parts[1]
      if (var_name %in% names(pf)) {
        if (length(parts) == 1) {
          # Case: "sex" (Just convert to factor)
          pf[[var_name]] <- as.factor(pf[[var_name]])
        } else if (length(parts) == 2) {
          # Case: "sex:male,female"
          levs <- unlist(strsplit(parts[2], ","))
          pf[[var_name]] <- factor(pf[[var_name]], levels = levs)
        } else if (length(parts) == 3) {
          # Case: "sex:1,2:male,female"
          levs   <- unlist(strsplit(parts[2], ","))
          labls  <- unlist(strsplit(parts[3], ","))
          pf[[var_name]] <- factor(pf[[var_name]], levels = levs, labels = labls)
        }
        if (verbose) { message(sprintf("Processed factor: %s", var_name))}
      } else {
        warning(sprintf("Variable '%s' not found in data frame.", var_name))
      }
    }
  }

  # Match DF data and NII data -------------------------------------------------
  ## 1. Parse id_var into a mapping (e.g., participant_id -> sub)
  if (verbose) { message("Matching participants in DF and NII") }
  id_map <- list()
  for (item in id_var) {
    parts <- unlist(strsplit(item, ":"))
    if (length(parts) == 2) {
      id_map[[parts[1]]] <- parts[2]
    } else {
      # Defaults for standard BIDS
      if (parts[1] == "participant_id") id_map[[parts[1]]] <- "sub"
      else if (parts[1] == "session_id") id_map[[parts[1]]] <- "ses"
      else id_map[[parts[1]]] <- parts[1] # Use var name as flag if not specified
    }
  }

  # 2. Get list of files from nii_data
  # nii_data can be a directory "path/to/data" or "path/to/data:_T1w"
  nii_parts <- unlist(strsplit(nii_data, ":"))
  search_dir <- nii_parts[1]
  file_suffix <- if (length(nii_parts) > 1) nii_parts[2] else ""
  if (dir.exists(search_dir)) {
    all_files <- list.files(search_dir, pattern = "\\.nii(\\.gz)?$",
                            full.names = TRUE, recursive = TRUE)
  } else {
    # If nii_data was already a list of files
    all_files <- nii_data
  }

  # 3. Match each row in 'pf' to exactly one file
  pf$nii_file <- as.character(NA)
  for (i in 1:nrow(pf)) {
    # Build a regex that matches ALL id_vars for this row
    # Example: sub-1234.*ses-5678.*aid-999.*_T1w.nii
    regex_parts <- c()
    for (var_name in names(id_map)) {
      flag <- id_map[[var_name]]
      val  <- pf[i, var_name]
      # BIDS flags are usually flag-value
      regex_parts <- c(regex_parts, sprintf("%s-%s", flag, val))
    }

    # Combine parts with .* to allow any flags in between
    # Add the user-specified suffix and the nifti extension
    match_pattern <- paste0(paste(regex_parts, collapse = ".*"),
                            ".*", file_suffix, "\\.nii(\\.gz)?$")

    matches <- all_files[grepl(match_pattern, all_files)]

    if (length(matches) == 0) {
      warning(sprintf("No matching NIfTI found for row %d (Pattern: %s)", i, match_pattern))
    } else if (length(matches) > 1) {
      stop(sprintf("Multiple files found for row %d: %s. Pattern was: %s",
                   i, paste(basename(matches), collapse=", "), match_pattern))
    } else {
      pf$nii_file[i] <- matches
    }
  }

  ## Remove rows where no NIfTI was found
  pf <- pf[!is.na(pf$nii_file), ]
  if (nrow(pf) == 0) { stop("Dataset is empty, please check inputs") }

  # Transfer & Decompress Participant Files ------------------------------------
  if (verbose) message(sprintf("Preparing %d files in scratch directory...", nrow(pf)))
  # We use prepNII to handle the heavy lifting (check existence, copy, gunzip)
  # but we use 'row_i' as the tag to ensure unique filenames in scratch.
  pf$nii_file <- sapply(1:nrow(pf), function(i) {
    # Call our standardized helper
    local_path <- prepNII(path = pf$nii_file[i],
                          tag = paste0("row_", i),
                          dir_scratch = dir_scratch)
    if (verbose && i %% 10 == 0) {
      message(sprintf("Prepared %d of %d files...", i, nrow(pf)))
    }
    return(local_path)
  })
  # Safety check for failed transfers
  if (any(is.na(pf$nii_file))) {
    stop("One or more participant NIfTI files could not be prepared in scratch.")
  }

  # load or generate mask -----------------------------------------------------
  if (!missing(roi_nii) && !is.null(roi_nii)) {
    if (verbose) message("Copying ROI mask to scratch")
    # Copy original file exactly as is to scratch
    local_roi <- file.path(dir_scratch, basename(roi_nii))
    file.copy(roi_nii, local_roi, overwrite = TRUE)
    # If it's zipped, decompress it locally on scratch
    if (grepl("\\.gz$", local_roi)) {
      roi_nii <- gunzip(local_roi, remove = TRUE, overwrite = TRUE)
    } else {
      roi_nii <- local_roi
    }
    # load mask
    if (verbose) message(sprintf("Loading ROI mask: %s", basename(roi_nii)))
    mask <- read.nii.volume(roi_nii,1)
    mask <- (mask != 0) * 1
    # gather image info
    img_dims <- info.nii(roi_nii, "dims")
    pixdim <- info.nii(roi_nii, "pixdim")
    orient <- info.nii(roi_nii, "orient")
  } else {
    if (verbose) message("No ROI specified. Generating mask from first data file...")
    mask <- read.nii.volume(pf$nii_file[1]) * 0 + 1
    # gather image info
    img_dims <- info.nii(pf$nii_file[1], "dims")
    pixdim   <- info.nii(pf$nii_file[1], "pixdim")
    orient   <- info.nii(pf$nii_file[1], "orient")
  }
  vxl_ls <- which(mask!=0, arr.ind=TRUE)

   # initialize log file if it doesn't exist --------------------------------------
  if (verbose) { message("Checking log file") }
  log.nii <- paste0(dir_scratch, "/log.nii")
  if (file.exists(log.nii) == FALSE || restart_log == TRUE) {
    if (verbose) { message("Continuing prior run") }
    init.nii(log.nii, dims=img_dims, pixdim=pixdim, orient=orient, init.value=0)
    write.nii.volume(log.nii, vol.num=1, value=mask)
  } else {
    if (verbose) { message("No log found, initializing...") }
    log <- read.nii.volume(log.nii,1)
    vxls_not_run <- (log == 1) * 1
    vxl_ls <- which(vxls_not_run==1, arr.ind=TRUE)
  }

  # set voxel looping parameters ------------------------------------------------
  n.vxls <- nrow(vxl_ls)
  ## check if there are no voxels
  if (n.vxls == 0) { stop("There are no voxels in the specified ROI to run") }
  if (verbose) { message(sprintf("Will process %d voxels.", n.vxls)) }

  ## randomize order ---
  if (rand_order) {
    if (verbose) { message("Randomizing voxel order") }
    vxl_ls <- vxl_ls[sample(1:n.vxls, n.vxls, replace=F), ]
  }

  # Check that voxels are valid --------------------------------------------------
  if (verbose) { message("Checking for invalid voxels") }
  valid_rows <- (vxl_ls[, 1] <= img_dims[1]) &
                (vxl_ls[, 2] <= img_dims[2]) &
                (vxl_ls[, 3] <= img_dims[3])
  valid_rows <- valid_rows & (vxl_ls[, 1] > 0) & (vxl_ls[, 2] > 0) & (vxl_ls[, 3] > 0)
  n_dropped <- sum(!valid_rows)
  if (n_dropped > 0) {
    warning(sprintf("Dropped %d voxels that were outside image boundaries (%d, %d, %d).",
                    n_dropped, img_dims[1], img_dims[2], img_dims[3]))
    vxl_ls <- vxl_ls[valid_rows, , drop = FALSE]
  }
  n.vxls <- nrow(vxl_ls)

  # specify model function -------------------------------------------------------
  ## -load voxelwise data
  ## -run user code (model_fcn)
  ## -save voxelwise output table
  ## -write log nii
  model.fxn <- function(coords, df, ...) {
    for (i in 1:nrow(df)) { df$nii[i] <- read.nii.voxel(df$nii_file[i], coords) }
    modelResult <- model_fcn(df)
    return(modelResult)
  }

  # run parallel loop over all voxels ----------------------------------------------
  do_debug=FALSE
  if (!is.na(debug) && debug > 0) {
      do_debug=TRUE
      message(sprintf("DEBUG MODE: Running first %d voxels sequentially...", debug))
      for (X in 1:debug) { model.fxn(X) }
      message("DEBUG DONE")
  } else {
    # Run voxels in parallel
    if (verbose) { message(sprintf("Starting voxelwise models on %d cores...", num_cores))}
    proc_start <- Sys.time()
    
    # Check Log for incomplete voxels ---------------------------------------------
    # initialize log file if it doesn't exist
    if (verbose) { message("Checking log file") }
    log.nii <- paste0(dir_scratch, "/log.nii")
    if (file.exists(log.nii) == FALSE) {
      if (verbose) { message("No log found, initializing...") }
      init.nii(log.nii, dims=img_dims, pixdim=pixdim, orient=orient, init.value=0)
      write.nii.volume(log.nii, vol.num=1, value=mask)
    } else if (restart_log == TRUE) {
      if (verbose) { message("Reinitializing log and restarting analysis...") }
      init.nii(log.nii, dims=img_dims, pixdim=pixdim, orient=orient, init.value=0)
      write.nii.volume(log.nii, vol.num=1, value=mask)
    } else {      
      if (verbose) { message("Continuing prior run") }
    }

    # read log file to get voxel list
    log <- read.nii.volume(log.nii,1)
    vxls_not_run <- (log == 1) * 1
    vxl_ls <- which(vxls_not_run==1, arr.ind=TRUE)
    n.vxls <- nrow(vxl_ls)
    if (n.vxls == 0) { stop("There are no voxels in the specified ROI to run") }
    if (verbose) { message(sprintf("Will process %d voxels.", n.vxls)) }
    # Check that voxels are valid --------------------------------------------------
    if (verbose) { message("Checking for invalid voxels") }
    valid_rows <- (vxl_ls[, 1] <= img_dims[1]) &
      (vxl_ls[, 2] <= img_dims[2]) &
      (vxl_ls[, 3] <= img_dims[3])
    valid_rows <- valid_rows & (vxl_ls[, 1] > 0) & (vxl_ls[, 2] > 0) & (vxl_ls[, 3] > 0)
    n_dropped <- sum(!valid_rows)
    if (n_dropped > 0) {
      warning(sprintf("Dropped %d voxels that were outside image boundaries (%d, %d, %d).",
                      n_dropped, img_dims[1], img_dims[2], img_dims[3]))
      vxl_ls <- vxl_ls[valid_rows, , drop = FALSE]
    }
    n.vxls <- nrow(vxl_ls)
    
    # repeat loop until all voxels are run -------------------------------------
    while (n.vxls > 0) {
      ## randomize voxel order ---
      if (rand_order) {
        if (verbose) { message("Randomizing voxel order") }
        vxl_ls <- vxl_ls[sample(1:n.vxls, n.vxls, replace=F), ]
      }
    
    # Split indices into a list of chunks
    chunker <- round(seq(1, n.vxls, length.out=num_cores+1))
    chunk_start <- chunker[1:(length(chunker)-1)]
    chunk_stop <- c(chunker[2:(length(chunker)-1)] - 1, chunker[length(chunker)])
    if (verbose) {
      for (i in 1:num_cores) {
        print(sprintf("worker_%02d: chunk %d to %d", i, chunk_start[i], chunk_stop[i]))
      }
    }
    registerDoParallel(num_cores)
    foreach(chk_id = 1:num_cores, .packages = all_libs, .errorhandling="pass") %dopar% {
      worker_start <- Sys.time()
      vxl_chunk <- vxl_ls[chunk_start[chk_id]:chunk_stop[chk_id], ]
      n_in_chunk <- nrow(vxl_chunk)
      worker_id <- sprintf("worker_%02d", chk_id)
      print(sprintf("starting worker_%02d", chk_id))
      for (i in 1:n_in_chunk) {
        tryCatch({
          coords <- vxl_chunk[i, ]
          df <- pf
          df$nii <- numeric(nrow(df))
          modelResult <- model.fxn(coords, df)
          ## add short delay periodically to reduce the likelihood of I/O collisions
          Sys.sleep(runif(1, 0, 0.01))
          table.to.nii(in.table=modelResult, coords=coords, save.dir=dir_scratch,
                       do.log=TRUE, model.string=model_pfx,
                       img.dims=img_dims, pixdim=pixdim, orient=orient)
          Sys.sleep(runif(1, 0, 0.01))
          write.nii.voxel(log.nii, coords, 2)
        }, error = function(e) {
          error_msg <- sprintf("Voxel %d failed: %s", X, e$message)
          write(error_msg, file = file.path(dir_scratch, "failed_voxels.log"), append = TRUE)
          write.nii.voxel(log.nii, coords, 3)
        })
        if (verbose) {
          pct <- floor((i / n_in_chunk) * 100)
          prev_pct <- floor(((i - 1) / n_in_chunk) * 100)
          if (pct > prev_pct) {
            message(sprintf("[%s] progress: %d%% completed", worker_id, pct))
          }
        }
      }
      worker_end <- Sys.time()
      elapsed <- difftime(worker_end, worker_start, units = "auto")
      message(sprintf("[%s] COMPLETED. Total elapsed time: %s", worker_id, format(elapsed)))
    }
    stopImplicitCluster() # Stop parallelization
    # read log file to get voxel list
    message("Checking log file to ensure completeness, will rerun missing voxels.")
    log <- read.nii.volume(log.nii,1)
    vxls_not_run <- (log == 1) * 1
    vxl_ls <- which(vxls_not_run==1, arr.ind=TRUE)
    n.vxls <- nrow(vxl_ls)
  }
  proc_stop <- Sys.time()
  elapsed <- difftime(proc_stop, proc_start, units = "auto")
  message(sprintf("MODELLING COMPLETED. Total elapsed time: %s", format(elapsed)))
  
  # Create Final Output Directory ------------------------------------------
  # We use the user's model prefix to create a specific sub-folder
  if (verbose) { message("Creating output directory") }
  output_dir <- file.path(dir_save, model_pfx)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Transfer Results from Scratch to Save ----------------------------------
  # Identify the results (everything in scratch except the temporary unzipped participant files)
  # -nifti.io will create files starting with the modelResult
  # -Compress in scratch first (faster local I/O).
  if (verbose) { message("Moving results and log files to final directory...") }
  result_files <- list.files(dir_scratch, pattern="modelResult", full.names=TRUE)
  for (f in result_files) {
    gz_file <- gzip(f, remove = TRUE, overwrite = TRUE)
    file.copy(gz_file, file.path(output_dir, basename(gz_file)), overwrite = TRUE)
  }
  gzip(file.path(dir_scratch, "log.nii"), remove = TRUE, overwrite = TRUE)
  file.copy(file.path(dir_scratch, "log.nii.gz"), output_dir, overwrite=TRUE)
  if (file.exists(file.path(dir_scratch, "failed_voxels.log"))) {
    file.copy(file.path(dir_scratch, "failed_voxels.log"), output_dir, overwrite=TRUE)
  }

  # Cleanup ----------------------------------------------------------------
  ## This removes the scratch folder and all unzipped participant copies
  if (cleanup) {
    if (verbose) message("Cleaning up scratch directory...")
    unlink(dir_scratch, recursive = TRUE)
  }
  if (verbose) message("Voxelwise analysis complete!")
  # Return the path to the results so the user can easily find them
  return(output_dir)
}
