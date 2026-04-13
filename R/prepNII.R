prepNII <- function(path, tag, dir_scratch) {
  # 1. Basic Safety Checks
  if (is.null(path) || length(path) == 0) return(NULL)

  # Ensure we only handle one file at a time
  if (length(path) > 1) {
    stop(sprintf("prepNII received a vector for %s. Expected a single path.", tag))
  }

  # Handle "none" string often used in masks
  if (path == "none") return(NULL)

  # 2. Define Destination
  # Ensure scratch directory exists
  if (!dir.exists(dir_scratch)) dir.create(dir_scratch, recursive = TRUE)

  # Strip .gz for the destination filename
  clean_name <- gsub("\\.gz$", "", basename(path))
  dest <- file.path(dir_scratch, paste0(tag, "_", clean_name))

  # 3. Process the File
  if (!file.exists(dest)) {
    if (grepl("\\.gz$", path)) {
      # Use remove = FALSE to keep the original .gz file in its source location
      R.utils::gunzip(path, destname = dest, remove = FALSE, overwrite = TRUE)
    } else {
      # If it's already .nii, just copy it to scratch
      file.copy(path, dest, overwrite = TRUE)
    }
  }

  return(dest)
}
