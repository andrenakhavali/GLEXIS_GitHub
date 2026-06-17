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
output_dir <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/Thresholds_wind/"
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
  
  
  
  wind_files <- tryCatch(get_var_files(hist_dir, "sfcwind"),
                         error = function(e) {
                           message(e$message)
                           stop("Failed to get sfcwind files")
                         })
  
  
  # Extract periods from filenames
  get_periods <- function(files) {
    periods <- sub(".*_(\\d{4}_\\d{4})\\.nc$", "\\1", files)
    if (any(is.na(periods))) {
      stop("Could not extract periods from filenames")
    }
    periods
  }
  periods_wind <- get_periods(wind_files)
  
  # Find common periods
  common_periods <- Reduce(intersect, list( periods_wind))
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
  
  wind_files <- filter_files(wind_files, common_periods)
  
  # Verify files exist
  check_files <- function(files) {
    for (f in files) {
      if (!file.exists(f)) {
        stop(paste("File does not exist:", f))
      }
    }
  }
  
  check_files(wind_files)
  
  # 2. Get Dimensions ------------------------------------------------------
  cat("\n===== Getting file dimensions =====\n")
  nc <- nc_open(wind_files[1])
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
  
  for (wind_file in wind_files) {
    nc <- nc_open(wind_file)
    time <- ncvar_get(nc, "time")
    time_origin <- as.Date(gsub("days since (.+)", "\\1", time_units))
    dates <- as.Date(time, origin = time_origin)
    years <- as.integer(format(dates, "%Y"))
    
    # Get days within the reference period
    ref_idx <- which(years >= ref_start & years <= ref_end)
    
    if (length(ref_idx) > 0) {
      ref_days <- ref_days + length(ref_idx)
      
      # Store reference indices for each file
      period <- sub(".*_(\\d{4}_\\d{4})\\.nc$", "\\1", wind_file)
      ref_days_list[[period]] <- list(
        file = wind_file,
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
    cat("- wind:", paste(periods_wind, collapse=", "), "\n")
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
  wind_ref <- array(NA_real_, dim = c(nlon, nlat, ref_days))
  
  # 5. Load Reference Data --------------------------------------------------
  cat("\n===== Loading reference data =====\n")
  current_idx <- 1
  
  for (period in names(ref_days_list)) {
    file_info <- ref_days_list[[period]]
    cat("\nProcessing period:", period, "\n")
    
    # Get corresponding files for this period
    wind_file <- grep(period, wind_files, value = TRUE)
    
    if (length(wind_file) != 1) {
      stop("Missing files for period ", period)
    }
    
    # Open files
    nc_wind <- nc_open(wind_file)
    
    # Calculate chunk size
    
    chunk_size <- file_info$count
    end_idx <- current_idx + chunk_size - 1
    
    # Read data in chunks with explicit numeric conversion
    cat("Reading wind...\n")
    wind_ref[,,current_idx:end_idx] <- as.numeric(ncvar_get(nc_wind, "sfcwind", 
                                                            start = c(1, 1, file_info$start_idx), 
                                                            count = c(-1, -1, file_info$count))) 
    
    
    # Close files
    nc_close(nc_wind)
    
    current_idx <- end_idx + 1
  }
  
  
  # 7. Calculate Variable-Specific Quantiles ---------------------------------
  cat("\n===== Calculating thresholds =====\n")
  
  quantile_sets <- list(
    wind = c(0.50, 0.75, 0.90, 0.95, 0.99))      # Upper tail for wind
  
  
  thresholds <- list(
    wind = array(NA_real_, dim = c(nlon, nlat, length(quantile_sets$wind)))
  )
  
  # Vectorized quantile calculation
  for (i in 1:nlon) {
    for (j in 1:nlat) {
      if (sum(!is.na(wind_ref[i,j,])) > 30) {
        thresholds$wind[i,j,] <- quantile(wind_ref[i,j,], quantile_sets$wind, na.rm = TRUE)
      }
    }
  }
  
  # 8. Save Results ---------------------------------------------------------
  output_file <- file.path(output_dir, 
                           paste0(model_pattern, "_thresholds_", 
                                  ref_start, "-", ref_end, ".nc"))
  
  # Create dimensions
  dim_lon <- ncdim_def("lon", "degrees_east", lon)
  dim_lat <- ncdim_def("lat", "degrees_north", lat)
  
  
  
  # Create output file name for wind only
  wind_output_file <- file.path(output_dir,
                                paste0(model_pattern, "_wind_thresholds_",
                                       ref_start, "-", ref_end, ".nc"))
  
  # Create dimensions
  dim_lon <- ncdim_def("lon", "degrees_east", lon)
  dim_lat <- ncdim_def("lat", "degrees_north", lat)
  dim_wind_perc <- ncdim_def("wind_percentile", "fraction", quantile_sets$wind)
  
  # Define only wind variables
  wind_var_defs <- list(
    sfcWind = ncvar_def("sfcWind", "m/s", list(dim_lon, dim_lat, dim_wind_perc),
                        longname = "Quantiles of daily surface wind speed",
                        prec = "float", missval = -9999),
    sfcWind_sd = ncvar_def("sfcWind_sd", "m/s", list(dim_lon, dim_lat),
                           longname = "Standard deviation of daily surface wind speed",
                           prec = "float", missval = -9999)
  )
  
  # Create NetCDF file with only wind data
  nc_wind <- nc_create(wind_output_file, wind_var_defs)
  
  safe_write <- function(data) {
    data[is.infinite(data)] <- -9999     # Replace infinite values with missing value code
    data[is.na(data)] <- -9999           # Replace NA values with missing value code
    return(data)
  }
  # Write wind data
  ncvar_put(nc_wind, "sfcWind", safe_write(thresholds$wind))
  
  # Add global attributes
  ncatt_put(nc_wind, 0, "title", paste("Wind Thresholds", model, ref_start, "-", ref_end))
  ncatt_put(nc_wind, 0, "percentiles", paste(quantile_sets$wind, collapse=", "))
  ncatt_put(nc_wind, 0, "missing_value", -9999)
  ncatt_put(nc_wind, 0, "created", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  nc_close(nc_wind)
  
  
  
  cat("\nProcessing complete for model", model, "! Results saved to:", wind_output_file, "\n")
}

cat("\n===== All models processed successfully! =====\n")

