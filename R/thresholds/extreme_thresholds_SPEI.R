# Load required libraries
library(ncdf4)
library(terra)
library(matrixStats)
library(SPEI)

# Configuration - MODIFY THESE
models <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")



ref_start <- 1990
ref_end <- 2010
min_year <- 1951
spei_scale <- 3

# Paths - UPDATE THESE
input_root <- "//pdrive/share/link/nakhavali.pdrv/watxene/ISIMIP/ISIMIP3b/InputData/climate_updated/bias-adjusted"
pet_root <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/PET"
output_dir <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/Thresholds_back/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Main processing loop for each model
for (model in models) {
  cat("\n\n===== Processing model:", model, "=====\n")
  
  # Model file pattern
  model_pattern <- switch(model,
                          "GFDL-ESM4" = "gfdl-esm4_r1i1p1f1_w5e5",
                          "IPSL-CM6A-LR" = "ipsl-cm6a-lr_r1i1p1f1_w5e5",
                          "MPI-ESM1-2-HR" = "mpi-esm1-2-hr_r1i1p1f1_w5e5",
                          "MRI-ESM2-0" = "mri-esm2-0_r1i1p1f1_w5e5",
                          "UKESM1-0-LL" = "ukesm1-0-ll_r1i1p1f2_w5e5"
  )
  
  # 1. Setup Data Sources ----------------------------------------------------
  cat("\n===== Setting up data sources =====\n")
  hist_dir <- file.path(input_root, "historical", model)
  pet_dir <- file.path(pet_root, "historical", model)
  
  # Robust file filtering function
  get_var_files <- function(dir, var) {
    all_files <- list.files(dir, 
                            pattern = paste0(model_pattern, ".*_", var, "_global_daily_.*\\.nc$"),
                            full.names = TRUE)
    
    if (length(all_files) == 0) {
      stop(paste("No", var, "files found in", dir))
    }
    
    # Extract years from filenames
    start_years <- as.integer(sub(".*_(\\d{4})_\\d{4}\\.nc$", "\\1", all_files))
    end_years <- as.integer(sub(".*_\\d{4}_(\\d{4})\\.nc$", "\\1", all_files))
    
    # Filter files that overlap with our reference period
    keep_files <- all_files[end_years >= ref_start & start_years <= ref_end]
    
    if (length(keep_files) == 0) {
      stop(paste("No", var, "files overlap with reference period", ref_start, "-", ref_end))
    }
    
    keep_files
  }
  
  # Get files for each variable with error handling
  tasmax_files <- tryCatch(get_var_files(hist_dir, "tasmax"),
                           error = function(e) {
                             message(e$message)
                             stop("Failed to get tasmax files")
                           })
  
  pr_files <- tryCatch(get_var_files(hist_dir, "pr"),
                       error = function(e) {
                         message(e$message)
                         stop("Failed to get pr files")
                       })
  
  tasmin_files <- tryCatch(get_var_files(hist_dir, "tasmin"),
                           error = function(e) {
                             message(e$message)
                             stop("Failed to get tasmin files")
                           })
  
  pet_files <- tryCatch(get_var_files(pet_dir, "pet"),
                        error = function(e) {
                          message(e$message)
                          stop("Failed to get pet files")
                        })
  
  # Verify we have files for all variables
  if (length(tasmax_files) == 0 || length(pr_files) == 0 || 
      length(tasmin_files) == 0 || length(pet_files) == 0) {
    cat("Skipping model", model, "due to missing files\n")
    next  # Skip to next model if files are missing
  }
  
  # Extract periods from filenames
  get_periods <- function(files) {
    periods <- sub(".*_(\\d{4}_\\d{4})\\.nc$", "\\1", files)
    if (any(is.na(periods))) {
      stop("Could not extract periods from filenames")
    }
    periods
  }
  
  periods_tasmax <- get_periods(tasmax_files)
  periods_pr <- get_periods(pr_files)
  periods_tasmin <- get_periods(tasmin_files)
  periods_pet <- get_periods(pet_files)
  
  # Find common periods
  common_periods <- Reduce(intersect, list(periods_tasmax, periods_pr, periods_tasmin, periods_pet))
  if (length(common_periods) == 0) {
    cat("No common periods for model", model, "- skipping\n")
    next
  }
  
  # Filter files to only include common periods
  filter_files <- function(files, periods) {
    result <- character(length(common_periods))
    for (i in seq_along(common_periods)) {
      match <- grep(common_periods[i], files, value = TRUE)
      if (length(match) == 0) {
        stop(paste("Missing file for period", common_periods[i]))
      }
      result[i] <- match[1]
    }
    result
  }
  
  tasmax_files <- filter_files(tasmax_files, common_periods)
  pr_files <- filter_files(pr_files, common_periods)
  tasmin_files <- filter_files(tasmin_files, common_periods)
  pet_files <- filter_files(pet_files, common_periods)
  
  # Verify files exist
  check_files <- function(files) {
    for (f in files) {
      if (!file.exists(f)) {
        stop(paste("File does not exist:", f))
      }
    }
  }
  
  check_files(tasmax_files)
  check_files(pr_files)
  check_files(tasmin_files)
  check_files(pet_files)
  
  # 2. Get Dimensions ------------------------------------------------------
  cat("\n===== Getting file dimensions =====\n")
  nc <- nc_open(tasmax_files[1])
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  nlon <- length(lon)
  nlat <- length(lat)
  time_units <- ncatt_get(nc, "time")$units
  nc_close(nc)
  
  cat("Grid dimensions:", nlon, "longitudes x", nlat, "latitudes\n")
  
  # 3. Calculate Reference Days --------------------------------------------
  cat("\n===== Calculating reference days =====\n")
  ref_days <- 0
  ref_days_list <- list()
  
  for (tasmax_file in tasmax_files) {
    nc <- nc_open(tasmax_file)
    time <- ncvar_get(nc, "time")
    time_origin <- as.Date(gsub("days since (.+)", "\\1", time_units))
    dates <- as.Date(time, origin = time_origin)
    years <- as.integer(format(dates, "%Y"))
    
    # Get days within the reference period
    ref_idx <- which(years >= ref_start & years <= ref_end)
    
    if (length(ref_idx) > 0) {
      ref_days <- ref_days + length(ref_idx)
      
      # Store reference indices for each file
      period <- sub(".*_(\\d{4}_\\d{4})\\.nc$", "\\1", tasmax_file)
      ref_days_list[[period]] <- list(
        file = tasmax_file,
        start_idx = min(ref_idx),
        count = length(ref_idx)
      )
    }
    nc_close(nc)
  }
  
  cat("Total reference days:", ref_days, "\n")
  
  if (ref_days == 0) {
    cat("\n===== Diagnostic Information =====\n")
    cat("Available periods in files:\n")
    cat("- tasmax:", paste(periods_tasmax, collapse=", "), "\n")
    cat("- pr:", paste(periods_pr, collapse=", "), "\n")
    cat("- tasmin:", paste(periods_tasmin, collapse=", "), "\n")
    cat("- pet:", paste(periods_pet, collapse=", "), "\n")
    
    cat("\nCommon periods:", paste(common_periods, collapse=", "), "\n")
    
    # Check specific years in files
    cat("\nChecking years in files:\n")
    for (file in tasmax_files[1:min(3, length(tasmax_files))]) {
      nc <- nc_open(file)
      time <- ncvar_get(nc, "time")
      time_origin <- as.Date(gsub("days since (.+)", "\\1", time_units))
      dates <- as.Date(time, origin = time_origin)
      years <- unique(as.integer(format(dates, "%Y")))
      nc_close(nc)
      cat("Years in", basename(file), ":", paste(range(years), collapse="-"), "\n")
    }
    
    cat("Skipping model", model, "due to no reference days\n")
    next
  }
  
  # 4. Initialize Arrays ----------------------------------------------------
  cat("\n===== Initializing arrays =====\n")
  tasmax_ref <- array(NA_real_, dim = c(nlon, nlat, ref_days))
  pr_ref <- array(NA_real_, dim = c(nlon, nlat, ref_days))
  tasmin_ref <- array(NA_real_, dim = c(nlon, nlat, ref_days))
  pet_ref <- array(NA_real_, dim = c(nlon, nlat, ref_days))
  
  # 5. Load Reference Data --------------------------------------------------
  cat("\n===== Loading reference data =====\n")
  current_idx <- 1
  
  for (period in names(ref_days_list)) {
    file_info <- ref_days_list[[period]]
    cat("\nProcessing period:", period, "\n")
    
    # Get corresponding files for this period
    tasmax_file <- grep(period, tasmax_files, value = TRUE)
    pr_file <- grep(period, pr_files, value = TRUE)
    tasmin_file <- grep(period, tasmin_files, value = TRUE)
    pet_file <- grep(period, pet_files, value = TRUE)
    
    if (length(tasmax_file) != 1 || length(pr_file) != 1 || 
        length(tasmin_file) != 1 || length(pet_file) != 1) {
      stop("Missing files for period ", period)
    }
    
    # Open files
    nc_tasmax <- nc_open(tasmax_file)
    nc_pr <- nc_open(pr_file)
    nc_tasmin <- nc_open(tasmin_file)
    nc_pet <- nc_open(pet_file)
    
    # Calculate chunk size
    chunk_size <- file_info$count
    end_idx <- current_idx + chunk_size - 1
    
    # Read data in chunks with explicit numeric conversion
    cat("Reading tasmax...\n")
    tasmax_ref[,,current_idx:end_idx] <- as.numeric(ncvar_get(nc_tasmax, "tasmax", 
                                                              start = c(1, 1, file_info$start_idx), 
                                                              count = c(-1, -1, file_info$count))) - 273.15
    
    cat("Reading pr...\n")
    pr_ref[,,current_idx:end_idx] <- as.numeric(ncvar_get(nc_pr, "pr", 
                                                          start = c(1, 1, file_info$start_idx), 
                                                          count = c(-1, -1, file_info$count))) * 86400
    
    cat("Reading tasmin...\n")
    tasmin_ref[,,current_idx:end_idx] <- as.numeric(ncvar_get(nc_tasmin, "tasmin", 
                                                              start = c(1, 1, file_info$start_idx), 
                                                              count = c(-1, -1, file_info$count))) - 273.15
    
    cat("Reading pet...\n")
    pet_ref[,,current_idx:end_idx] <- as.numeric(ncvar_get(nc_pet, "pet", 
                                                           start = c(1, 1, file_info$start_idx), 
                                                           count = c(-1, -1, file_info$count)))
    
    # Close files
    nc_close(nc_tasmax)
    nc_close(nc_pr)
    nc_close(nc_tasmin)
    nc_close(nc_pet)
    
    current_idx <- end_idx + 1
  }
  
  # 6. SPEI Calculation (Optimized) ----------------------------------------
  cat("\n===== Calculating SPEI-", spei_scale, " =====\n", sep="")
  
  # Convert to 2D matrix (grid cells x time)
  wb_daily <- pr_ref - pet_ref
  wb_matrix <- matrix(wb_daily, nrow = nlon*nlat, ncol = dim(wb_daily)[3])
  spei_matrix <- matrix(NA_real_, nrow = nlon*nlat, ncol = dim(wb_daily)[3])
  
  # Progress tracking
  cat("Processing", nlon*nlat, "grid cells...\n")
  start_time <- Sys.time()
  
  for (i in 1:nrow(wb_matrix)) {
    # Progress update every 5%
    if (i %% floor(nrow(wb_matrix)/20) == 0) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      cat(sprintf("Progress: %d%%, Elapsed: %.1f mins\n",
                  round(100*i/nrow(wb_matrix)), elapsed))
    }
    
    if (sum(!is.na(wb_matrix[i,])) > 30) {  # Minimum 30 days of data
      suppressWarnings({
        spei_matrix[i,] <- spei(wb_matrix[i,],
                                scale = spei_scale,
                                kernel = list(type = "rectangular", shift = 0),
                                verbose = FALSE)$fitted
      })
    }
  }
  
  # Convert back to 3D array
  spei_values <- array(spei_matrix, dim = dim(wb_daily))
  
  # 7. Calculate Variable-Specific Quantiles ---------------------------------
  cat("\n===== Calculating thresholds =====\n")
  
  quantile_sets <- list(
    tasmax = c(0.50, 0.75, 0.90, 0.95, 0.99),  # Upper tail for heat
    pr = c(0.50, 0.75, 0.90, 0.95, 0.99),      # Upper tail for rain
    tasmin = c(0.01, 0.05, 0.10, 0.25, 0.50),  # Lower tail for cold
    spei = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)  # Full range
  )
  
  thresholds <- list(
    tasmax = array(NA_real_, dim = c(nlon, nlat, length(quantile_sets$tasmax))),
    pr = array(NA_real_, dim = c(nlon, nlat, length(quantile_sets$pr))),
    tasmin = array(NA_real_, dim = c(nlon, nlat, length(quantile_sets$tasmin))),
    spei = array(NA_real_, dim = c(nlon, nlat, length(quantile_sets$spei)))
  )
  
  # Vectorized quantile calculation
  for (i in 1:nlon) {
    for (j in 1:nlat) {
      if (sum(!is.na(tasmax_ref[i,j,])) > 30) {
        thresholds$tasmax[i,j,] <- quantile(tasmax_ref[i,j,], quantile_sets$tasmax, na.rm = TRUE)
      }
      if (sum(!is.na(pr_ref[i,j,])) > 30) {
        thresholds$pr[i,j,] <- quantile(pr_ref[i,j,], quantile_sets$pr, na.rm = TRUE)
      }
      if (sum(!is.na(tasmin_ref[i,j,])) > 30) {
        thresholds$tasmin[i,j,] <- quantile(tasmin_ref[i,j,], quantile_sets$tasmin, na.rm = TRUE)
      }
      if (sum(!is.na(spei_values[i,j,])) > 30) {
        thresholds$spei[i,j,] <- quantile(spei_values[i,j,], quantile_sets$spei, na.rm = TRUE)
      }
    }
  }
  
  # 8. Save Results ---------------------------------------------------------
  output_file <- file.path(output_dir, 
                           paste0(model_pattern, "_thresholds_", 
                                  ref_start, "-", ref_end, "_SPEI", spei_scale, ".nc"))
  
  # Create dimensions
  dim_lon <- ncdim_def("lon", "degrees_east", lon)
  dim_lat <- ncdim_def("lat", "degrees_north", lat)
  
  # Variable-specific percentile dimensions
  dim_tasmax_perc <- ncdim_def("tasmax_percentile", "fraction", quantile_sets$tasmax)
  dim_pr_perc <- ncdim_def("pr_percentile", "fraction", quantile_sets$pr)
  dim_tasmin_perc <- ncdim_def("tasmin_percentile", "fraction", quantile_sets$tasmin)
  dim_spei_perc <- ncdim_def("spei_percentile", "fraction", quantile_sets$spei)
  
  # Define variables with explicit range limits
  var_defs <- list(
    tasmax = ncvar_def("tasmax", "degC", list(dim_lon, dim_lat, dim_tasmax_perc),
                       longname = "Quantiles of daily maximum temperature",
                       prec = "float", missval = -9999),
    pr = ncvar_def("pr", "mm/day", list(dim_lon, dim_lat, dim_pr_perc),
                   longname = "Quantiles of daily precipitation",
                   prec = "float", missval = -9999),
    tasmin = ncvar_def("tasmin", "degC", list(dim_lon, dim_lat, dim_tasmin_perc),
                       longname = "Quantiles of daily minimum temperature",
                       prec = "float", missval = -9999),
    spei = ncvar_def("spei", "unitless", list(dim_lon, dim_lat, dim_spei_perc),
                     longname = paste("SPEI-", spei_scale, " quantiles"),
                     prec = "float", missval = -9999)
  )
  
  # Create NetCDF file
  nc <- nc_create(output_file, var_defs)
  
  # Handle potential infinite/NA values before writing
  safe_write <- function(data) {
    data[is.infinite(data)] <- -9999
    data[is.na(data)] <- -9999
    return(data)
  }
  
  # Write data with safety checks
  ncvar_put(nc, "tasmax", safe_write(thresholds$tasmax))
  ncvar_put(nc, "pr", safe_write(thresholds$pr))
  ncvar_put(nc, "tasmin", safe_write(thresholds$tasmin))
  ncvar_put(nc, "spei", safe_write(thresholds$spei))
  
  # Add global attributes
  ncatt_put(nc, 0, "title", paste("Climate Thresholds", model, ref_start, "-", ref_end))
  ncatt_put(nc, 0, "SPEI_scale", spei_scale)
  ncatt_put(nc, 0, "missing_value", -9999)
  ncatt_put(nc, 0, "created", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  nc_close(nc)
  gc()
  cat("\nProcessing complete for model", model, "! Results saved to:", output_file, "\n")
}

cat("\n===== All models processed successfully! =====\n")

