# ======================================================================
# GLEXIS: Global Extreme Events from ISIMIP Simulations
# A comprehensive processing pipeline for ISIMIP3b climate extremes
# Developed by Dr. Andre Nakhavali (nakhavali@iiasa.ac.at)
# International Institute for Applied Systems Analysis (IIASA)
# License: Creative Commons Attribution 4.0 International (CC BY 4.0)
# ======================================================================

# Load required libraries
library(ncdf4)
library(terra)
library(SPEI)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lubridate)
library(yaml)
library(abind)
library(exactextractr)
library(tidyverse)
library(parallel)

# =============== CONFIGURATION ===============
# Package metadata
pkg_metadata <- list(
  name = "GLEXIS",
  version = "1.0.0",
  description = "Global Extreme Events from ISIMIP Simulations",
  author = "Dr. Andre Nakhavali",
  email = "nakhavali@iiasa.ac.at",
  institution = "International Institute for Applied Systems Analysis (IIASA)",
  license = "CC BY 4.0",
  doi = "10.5281/zenodo.20734236",
  reference = "Nakhavali et al. (2025) GLEXIS: A Global Extreme Events Dataset from ISIMIP3b"
)

# Model and scenario configuration
models <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
scenarios <- c("ssp126", "ssp370", "ssp585")
years <- 2015:2100

# Define the decades (20-year chunks ending at 2100)
decade_starts <- seq(2015, 2095, by = 20)  # 2015, 2035, 2055, 2075, 2095
decade_ends <- pmin(decade_starts + 19, 2100)  # Cap at 2100

# Create decade labels
decade_labels <- paste(decade_starts, decade_ends, sep = "-")

# Create all combinations
combinations <- expand.grid(
  model = models,
  scenario = scenarios,
  start = decade_starts,
  stringsAsFactors = FALSE
) %>%
  mutate(
    end = decade_ends[match(start, decade_starts)]
  ) %>%
  select(model, scenario, start, end) %>%
  arrange(model, scenario, start)

args <- commandArgs(trailingOnly=TRUE)
scenario_index <- as.integer(args[1])

# Create station mapping
station <- list()
x <- 0
for (i in 1:75) {
  if (i > 1) { x <- x + 1 }
  station[[i]] <- (1 + x)
}

################################

s <- 1#station[[scenario_index]]
model <- combinations$model[s]
scenario <- combinations$scenario[s]
years <- combinations$start[s]:combinations$end[s]

# Core processing parameters
params <- list(
  run_type = "main",  # "main" or "ensemble"
  spei_scale = 3,
  south_pole_cutoff = -60,
  min_precip_filter = 1,
  hyper_arid_threshold = 5,
  annual_arid_threshold = 30,
  annual_wet_threshold = 2000,
  wind_min_speed = 8,
  tropical_lat_limit = 23.5,
  equatorial_lat_limit = 5,
  jan_precip_threshold = 10,
  baseline_period = "1990-2010",
  
  # MAIN RUN THRESHOLDS
  thresholds_quantiles = list(
    tasmax = 4,      # 95th percentile
    pr = 4,          # 95th percentile
    tasmin = 1,      # p01 threshold-file layer
    sfcwind = 1,     # Only one threshold in separated files
    windgust = 1,    # Only one threshold in separated files
    spei = 2         # p05 baseline SPEI-3 threshold-file layer
  ),
  
  # ENSEMBLE RUN THRESHOLDS
  ensemble_thresholds = list(
    tasmax = c(5, 3, 2),  # 99th, 90th, 75th percentiles
    pr = c(5, 3, 2),
    tasmin = c(2, 3, 5),
    sfcwind = c(5, 3, 2),
    windgust = c(5, 3, 2),
    spei = c(3, 4, 5)
  ),
  
  targeted_fixes = list(
    atacama_rainfall = TRUE,
    caribbean_cold = TRUE,
    ethiopia_heat = FALSE,
    north_atlantic_wind = TRUE,
    caribbean_drought = TRUE
  )
)

# File chunk structure
file_chunks <- list(
  "2015_2020" = 2015:2020,
  "2021_2030" = 2021:2030,
  "2031_2040" = 2031:2040,
  "2041_2050" = 2041:2050,
  "2051_2060" = 2051:2060,
  "2061_2070" = 2061:2070,
  "2071_2080" = 2071:2080,
  "2081_2090" = 2081:2090,
  "2091_2100" = 2091:2100
)

# Path configuration
input_root <- "//pdrive/share/link/nakhavali.pdrv/watxene/ISIMIP/ISIMIP3b/InputData/climate_updated/bias-adjusted"
pet_root <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/PET"
thresholds_dir <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/Thresholds_back/"
wind_thresholds_dir <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/Thresholds_separated/"  
output_dir <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/GLEXIS_v2/"

# Create directory structure
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "annual"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "decadal"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "ensemble"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "regional"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "ancillary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "documentation"), recursive = TRUE, showWarnings = FALSE)

# =============== GLOBAL ATTRIBUTES ===============
global_attrs <- list(
  title = "GLEXIS: Global Extreme Events from ISIMIP Simulations",
  summary = paste("Monthly counts of extreme climate events (heat, cold, precipitation,",
                  "wind, drought) derived from ISIMIP3b bias-adjusted climate projections",
                  "for multiple models and scenarios."),
  keywords = paste("climate extremes, ISIMIP3b, heat waves, cold spells, heavy precipitation,",
                   "drought, wind storms, climate change impacts"),
  institution = pkg_metadata$institution,
  source = "ISIMIP3b (https://www.isimip.org/) and W5E5 (https://doi.org/10.5880/PIK.2019.023)",
  references = paste("Lange (2021) https://doi.org/10.5194/essd-13-2079-2021;",
                     pkg_metadata$reference),
  license = pkg_metadata$license,
  Conventions = "CF-1.10",
  comment = paste("Processed using GLEXIS R package with targeted fixes for known artifacts.",
                  "See documentation for exclusion mask details and methodology."),
  creator_name = pkg_metadata$author,
  creator_email = pkg_metadata$email,
  creator_institution = pkg_metadata$institution,
  date_created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
  software_name = pkg_metadata$name,
  software_version = pkg_metadata$version,
  doi = pkg_metadata$doi
)

# =============== SUPPORT FUNCTIONS ===============

# Vectorized wind gust estimation function
estimate_wind_gust_vectorized <- function(mean_wind, temperature, humidity, pressure, z0 = 0.03) {
  # Constants
  k <- 0.4   # von Karman constant
  z <- 10    # Measurement height (m)
  g <- 9.81  # Gravity (m/s²)
  Rd <- 287.05
  Rv <- 461.5
  Cp <- 1004
  Lv <- 2.5e6
  
  # Calculate saturation vapor pressure (Pa)
  es <- 6.112 * exp(17.67 * temperature / (temperature + 243.5)) * 100
  e <- humidity * es / 100
  rho <- pressure / (Rd * temperature) * (1 - 0.378 * e / pressure)
  
  # Friction velocity
  u_star <- mean_wind * k / log(z / z0)
  
  # Sensible heat flux approximation
  H <- 0.1 * rho * Cp * u_star * (temperature - 273.15)
  
  # Obukhov length
  L <- - (u_star^3 * rho * temperature) / (k * g * H)
  
  # Stability parameter
  zeta <- z / L
  
  # Gust factor based on stability
  gust_factor <- rep(1.45, length(zeta))
  
  # Vectorized stability conditions
  very_unstable <- zeta < -0.5
  mod_unstable <- zeta >= -0.5 & zeta < 0
  mod_stable <- zeta > 0 & zeta <= 0.5
  very_stable <- zeta > 0.5
  
  gust_factor[very_unstable] <- 1.6 + 0.1 * abs(zeta[very_unstable])^0.25
  gust_factor[mod_unstable] <- 1.5 + 0.1 * abs(zeta[mod_unstable])^0.3
  gust_factor[mod_stable] <- 1.4 - 0.15 * zeta[mod_stable]^0.35
  gust_factor[very_stable] <- 1.2 - 0.1 * zeta[very_stable]^0.4
  
  # Apply bounds
  gust_factor <- pmax(1.2, pmin(2.0, gust_factor))
  
  # Calculate wind gust
  wind_gust <- mean_wind * gust_factor
  return(wind_gust)
}

# Function to get model-specific file pattern
get_model_pattern <- function(model) {
  if (model == "UKESM1-0-LL") return("ukesm1-0-ll_r1i1p1f2_w5e5")
  switch(model,
         "GFDL-ESM4" = "gfdl-esm4_r1i1p1f1_w5e5",
         "IPSL-CM6A-LR" = "ipsl-cm6a-lr_r1i1p1f1_w5e5", 
         "MPI-ESM1-2-HR" = "mpi-esm1-2-hr_r1i1p1f1_w5e5",
         "MRI-ESM2-0" = "mri-esm2-0_r1i1p1f1_w5e5",
         "UKESM1-0-LL" = "ukesm1-0-ll_r1i1p1f2_w5e5"
  )
}

# Function to create CF-compliant dimensions
create_cf_dims <- function(lons, lats, times) {
  list(
    lon = ncdim_def("lon", "degrees_east", lons, longname = "longitude"),
    lat = ncdim_def("lat", "degrees_north", lats, longname = "latitude"),
    time = ncdim_def("time", "days since 1850-01-01", 
                     as.numeric(times - as.Date("1850-01-01")),
                     calendar = "proleptic_gregorian",
                     longname = "time")
  )
}

# Function to save ancillary outputs
save_ancillary_outputs <- function(model, scenario, year, ancillary_data, lons, lats, dates) {
  cat("\nSaving ancillary outputs for", model, scenario, year, "...\n")
  
  # 1. Create output directory with proper permissions
  anc_dir <- file.path(output_dir, "ancillary", paste0(tolower(model), "_", scenario), year)
  
  # First check if directory exists and is writable
  if (dir.exists(anc_dir)) {
    if (file.access(anc_dir, 2) != 0) { # 2 = write permission
      stop(paste("No write permission in existing directory:", anc_dir))
    }
  } else {
    # Try to create directory with different permission modes
    tryCatch({
      dir.create(anc_dir, recursive = TRUE, mode = "0775", showWarnings = FALSE)
    }, warning = function(w) {
      cat("Warning creating directory:", w$message, "\n")
    }, error = function(e) {
      stop(paste("Failed to create directory:", anc_dir, "\nError:", e$message))
    })
    
    # Verify directory was created
    if (!dir.exists(anc_dir)) {
      # Try alternative location if primary fails
      alt_dir <- file.path(tempdir(), "GLEXIS_ancillary", 
                           paste0(tolower(model), "_", scenario), year)
      dir.create(alt_dir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(alt_dir)) {
        stop("Failed to create directory in both primary and temp locations")
      }
      anc_dir <- alt_dir
      warning(paste("Using temporary directory:", anc_dir))
    }
  }
  
  # 2. Prepare NetCDF file with robust error handling
  nc_file <- file.path(anc_dir, 
                       paste0("GLEXIS_", tolower(model), "_", scenario, "_", year,
                              "_ancillary.nc"))
  
  # 3. Create dimensions and variables
  dims <- list(
    lon = ncdim_def("lon", "degrees_east", lons),
    lat = ncdim_def("lat", "degrees_north", lats),
    time = ncdim_def("time", "days since 1850-01-01", 
                     as.numeric(dates - as.Date("1850-01-01")),
                     calendar = "standard")
  )
  
  vars <- list()
  for (name in names(ancillary_data)) {
    # Ensure data is proper 3D array
    if (!is.array(ancillary_data[[name]]) || length(dim(ancillary_data[[name]])) != 3) {
      ancillary_data[[name]] <- array(ancillary_data[[name]], 
                                      dim = c(length(lons), length(lats), length(dates)))
    }
    
    vars[[name]] <- ncvar_def(
      name = name,
      units = "1",
      dim = dims,
      missval = -9999,
      longname = name,
      prec = "float",
      compression = 5,
      chunksizes = c(min(50, length(lons)), min(50, length(lats)), 1)
    )
  }
  
  # 4. Save data with multiple fallback options
  tryCatch({
    # Attempt 1: Regular save
    nc <- nc_create(nc_file, vars, force_v4 = TRUE)
    for (name in names(ancillary_data)) {
      ncvar_put(nc, name, aperm(ancillary_data[[name]], c(1, 2, 3)))
    }
    nc_close(nc)
    cat("Successfully saved ancillary outputs:", nc_file, "\n")
  }, error = function(e) {
    cat("Primary save failed, attempting fallback methods...\n")
    
    # Attempt 2: Save to temporary location
    temp_file <- file.path(tempdir(), basename(nc_file))
    tryCatch({
      nc <- nc_create(temp_file, vars, force_v4 = TRUE)
      for (name in names(ancillary_data)) {
        ncvar_put(nc, name, aperm(ancillary_data[[name]], c(1, 2, 3)))
      }
      nc_close(nc)
      
      # Try to move to final location
      if (file.copy(temp_file, nc_file)) {
        unlink(temp_file)
        cat("Saved via temp file copy:", nc_file, "\n")
      } else {
        warning(paste("Could not move temp file to final location. File remains at:", temp_file))
      }
    }, error = function(e2) {
      # Attempt 3: Save as RDS as last resort
      rds_file <- sub("\\.nc$", ".rds", nc_file)
      saveRDS(list(data = ancillary_data, dims = dims, vars = vars), rds_file)
      warning(paste("NetCDF save failed. Saved as RDS:", rds_file))
    })
  })
  
  # 5. Verify file was created
  if (!file.exists(nc_file) && !file.exists(sub("\\.nc$", ".rds", nc_file))) {
    stop("Failed to save ancillary outputs in any format")
  }
}

# Function to save regional outputs
save_regional_outputs <- function(model, scenario, year, monthly_event_counts, lons, lats) {
  cat("\nSaving regional outputs with validation...\n")
  
  # Create output directory
  reg_dir <- file.path(output_dir, "regional", paste0(tolower(model), "_", scenario), year)
  dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load country boundaries
  countries <- ne_countries(scale = 50, returnclass = "sf")
  
  # Create raster template matching climate grid
  r <- rast(
    ncols = length(lons),
    nrows = length(lats),
    xmin = min(lons),
    xmax = max(lons),
    ymin = min(lats),
    ymax = max(lats),
    crs = "EPSG:4326"
  )
  
  # Process each event type
  for (event in names(monthly_event_counts)) {
    # Validate input data first
    if (all(monthly_event_counts[[event]] == 0)) {
      warning(paste("All zero values detected for", event, "- checking thresholds"))
      
      # Diagnostic output
      cat("\nEvent counts validation for", event, ":\n")
      print(summary(as.vector(monthly_event_counts[[event]])))
    }
    
    # Create monthly rasters with proper orientation
    monthly_rasters <- lapply(1:12, function(m) {
      # Transpose to match raster orientation (longitude in columns)
      vals <- t(monthly_event_counts[[event]][, , m])
      r <- rast(vals)
      dx <- median(diff(sort(unique(lons))))
      dy <- median(diff(sort(unique(lats))))
      ext(r) <- ext(c(min(lons) - dx / 2, max(lons) + dx / 2,
                      min(lats) - dy / 2, max(lats) + dy / 2))
      crs(r) <- "EPSG:4326"
      return(r)
    })
    
    # Aggregate by country for each month
    results <- list()
    for (m in 1:12) {
      # Use exact_extract for precise area-weighted sums
      country_stats <- exact_extract(monthly_rasters[[m]], countries, 'sum')
      
      results[[m]] <- data.frame(
        ISO3 = countries$iso_a3,
        Country = countries$name,
        Event = event,
        Year = year,
        Month = m,
        Count = round(country_stats, 1)  # Round to 1 decimal place
      )
    }
    
    # Combine all months
    results_df <- do.call(rbind, results)
    
    # Verify we have non-zero counts
    if (all(results_df$Count == 0)) {
      warning(paste("All zero counts for", event, "in", year))
      
      # Save diagnostic information
      diag_file <- file.path(reg_dir, 
                             paste0("DIAGNOSTIC_", event, "_", year, ".txt"))
      cat("Diagnostic for", event, year, ":\n",
          "Min count:", min(results_df$Count), "\n",
          "Max count:", max(results_df$Count), "\n",
          "Mean count:", mean(results_df$Count), "\n",
          file = diag_file)
    }
    
    # Save to CSV
    csv_file <- file.path(reg_dir, 
                          paste0("GLEXIS_", tolower(model), "_", scenario, "_", year,
                                 "_", event, "_country_counts.csv"))
    write.csv(results_df, csv_file, row.names = FALSE)
    cat("Saved regional outputs:", csv_file, "\n")
    
    # Create diagnostic plot for first month
    png(file.path(reg_dir, paste0("diagnostic_", event, "_month1.png")))
    plot(monthly_rasters[[1]], main = paste(event, "Month 1"))
    plot(st_geometry(countries), add = TRUE, border = "red")
    dev.off()
  }
}

# =============== DECADAL AGGREGATION FUNCTION ===============
aggregate_decadal_outputs <- function() {
  cat("\n===== AGGREGATING DECADAL OUTPUTS =====\n")
  
  # Define decades
  decades <- list(
    "2021_2030" = 2021:2030,
    "2031_2040" = 2031:2040,
    "2041_2050" = 2041:2050,
    "2051_2060" = 2051:2060,
    "2061_2070" = 2061:2070,
    "2071_2080" = 2071:2080,
    "2081_2090" = 2081:2090,
    "2091_2100" = 2091:2100
  )
  
  # Process each model, scenario and decade
  for (model in models) {
    for (scenario in scenarios) {
      for (decade_name in names(decades)) {
        decade_years <- decades[[decade_name]]
        cat("\nAggregating", model, scenario, decade_name, "...\n")
        
        # Create output directory
        dec_dir <- file.path(output_dir, "decadal", paste0(tolower(model), "_", scenario), decade_name)
        dir.create(dec_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Initialize decade rasters
        decade_rasters <- list()
        
        # Check if we have all years
        all_years_present <- TRUE
        for (year in decade_years) {
          year_dir <- file.path(output_dir, "annual", paste0(tolower(model), "_", scenario), year)
          if (!dir.exists(year_dir)) {
            cat("Missing year:", year, "\n")
            all_years_present = FALSE
            break
          }
        }
        
        if (!all_years_present) {
          cat("Skipping decade", decade_name, "due to missing years\n")
          next
        }
        
        # Aggregate each event type
        for (event in c("heat", "rain", "cold", "wind", "drought")) {
          # Initialize decade array
          decade_array <- NULL
          
          for (year in decade_years) {
            year_dir <- file.path(output_dir, "annual", paste0(tolower(model), "_", scenario), year)
            tif_file <- file.path(year_dir, 
                                  paste0("GLEXIS_", tolower(model), "_", scenario, "_", year, 
                                         "_", event, "_monthly_counts.tif"))
            
            if (file.exists(tif_file)) {
              r <- rast(tif_file)
              if (is.null(decade_array)) {
                decade_array <- array(0, dim = c(dim(r), length(decade_years)))
              }
              decade_array[, , which(decade_years == year)] <- as.array(r)
            } else {
              cat("Missing file:", tif_file, "\n")
              all_years_present = FALSE
              break
            }
          }
          
          if (!all_years_present) break
          
          # Calculate decade average
          decade_avg <- apply(decade_array, c(1,2,3), mean, na.rm = TRUE)
          
          # Create raster
          r_decade <- rast(decade_avg)
          ext(r_decade) <- ext(r)
          crs(r_decade) <- crs(r)
          
          # Save decade average
          out_file <- file.path(dec_dir, 
                                paste0("GLEXIS_", tolower(model), "_", scenario, "_", decade_name, 
                                       "_", event, "_monthly_avg.tif"))
          
          # Create raster with metadata
          r_decade <- rast(decade_avg)
          ext(r_decade) <- ext(r)
          crs(r_decade) <- crs(r)
          names(r_decade) <- month.name[1:12]
          time(r_decade) <- as.Date(paste0(substr(decade_name, 1, 4), "-", 1:12, "-15"))
          varnames(r_decade) <- "event_days"
          units(r_decade) <- "days"
          
          # Set long names for each layer
          for (i in 1:12) {
            longnames(r_decade)[i] <- paste("Decadal average", event, "days in", month.name[i], decade_name)
          }
          
          writeRaster(r_decade, filename = out_file, overwrite = TRUE,
                      gdal = c("DESCRIPTION=Decadal average monthly event counts"))
          cat("Saved:", out_file, "\n")
          
          # Store for ensemble
          decade_rasters[[event]] <- r_decade
        }
        
        # Store decade results for ensemble processing
        if (all_years_present) {
          ensemble_data[[paste(model, scenario, decade_name)]] <- decade_rasters
        }
      }
    }
  }
}

# =============== ENSEMBLE PROCESSING FUNCTION ===============
process_ensemble_outputs <- function() {
  cat("\n===== PROCESSING ENSEMBLE OUTPUTS =====\n")
  
  # Define decades
  decades <- names(file_chunks)[-1]  
  
  # Process each scenario and decade
  for (scenario in scenarios) {
    for (decade_name in decades) {
      cat("\nProcessing ensemble for", scenario, decade_name, "...\n")
      
      # Create output directory
      ens_dir <- file.path(output_dir, "ensemble", scenario, decade_name)
      dir.create(ens_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Collect all models for this scenario and decade
      model_rasters <- list()
      for (model in models) {
        key <- paste(model, scenario, decade_name)
        if (key %in% names(ensemble_data)) {
          model_rasters[[model]] <- ensemble_data[[key]]
        }
      }
      
      if (length(model_rasters) == 0) {
        cat("No models found for", scenario, decade_name, "\n")
        next
      }
      
      # Process each event type
      for (event in c("heat", "rain", "cold", "wind", "drought")) {
        # Create stack of all models for this event
        event_stack <- list()
        for (model in names(model_rasters)) {
          if (event %in% names(model_rasters[[model]])) {
            event_stack[[model]] <- model_rasters[[model]][[event]]
          }
        }
        
        if (length(event_stack) == 0) next
        
        # Create raster stack
        r_stack <- rast(event_stack)
        
        # Calculate ensemble statistics
        ens_stats <- list(
          mean = app(r_stack, fun = mean, na.rm = TRUE),
          sd = app(r_stack, fun = sd, na.rm = TRUE),
          min = app(r_stack, fun = min, na.rm = TRUE),
          max = app(r_stack, fun = max, na.rm = TRUE),
          p10 = app(r_stack, fun = function(x) quantile(x, probs = 0.1, na.rm = TRUE)),
          p50 = app(r_stack, fun = function(x) quantile(x, probs = 0.5, na.rm = TRUE)),
          p90 = app(r_stack, fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE))
        )
        
        # Save each statistic
        for (stat_name in names(ens_stats)) {
          out_file <- file.path(ens_dir, 
                                paste0("GLEXIS_ensemble_", scenario, "_", decade_name, 
                                       "_", event, "_", stat_name, ".tif"))
          
          # Create raster with metadata
          r_stat <- ens_stats[[stat_name]]
          names(r_stat) <- month.name[1:12]
          time(r_stat) <- as.Date(paste0(substr(decade_name, 1, 4), "-", 1:12, "-15"))
          varnames(r_stat) <- "event_days"
          units(r_stat) <- "days"
          
          # Set long names for each layer
          for (i in 1:12) {
            longnames(r_stat)[i] <- paste("Ensemble", stat_name, "of", event, "days in", month.name[i], decade_name)
          }
          
          writeRaster(r_stat, filename = out_file, overwrite = TRUE,
                      gdal = c("DESCRIPTION=Ensemble statistics for monthly event counts"))
          cat("Saved:", out_file, "\n")
        }
      }
    }
  }
}

# =============== DOCUMENTATION GENERATION ===============
generate_documentation <- function() {
  cat("\n===== GENERATING DOCUMENTATION =====\n")
  
  # Create README
  readme_content <- paste(
    "GLEXIS: Global Extreme Events from ISIMIP Simulations\n",
    "=====================================================\n",
    "\n",
    "Dataset Description:\n",
    "--------------------\n",
    "This dataset contains monthly counts of extreme climate events derived from ISIMIP3b ",
    "bias-adjusted climate projections for multiple models and scenarios. The dataset ",
    "includes five daily threshold indicators: hot days, heavy-rain days, cold days, ",
    "high-wind days, and standardized water-balance lower-tail days.\n",
    "\n",
    "Methodology:\n",
    "-----------\n",
    "1. Extreme events are defined using percentile thresholds from a 1990-2010 baseline period\n",
    "2. Heat days: Daily maximum temperature > 95th percentile\n",
    "3. Rain days: Daily precipitation > 95th percentile\n",
    "4. Cold days: Daily minimum temperature < 1st percentile\n",
    "5. Wind days: Daily wind speed > 95th percentile (mean wind) or > 99th percentile (max wind)\n",
    "6. Water-balance lower-tail days: Annual within-cell standardized precipitation minus PET < the 5th-percentile threshold-file layer\n",
    "\n",
    "Data Structure:\n",
    "---------------\n",
    "The dataset is organized into five main directories:\n",
    "1. annual/: Contains annual files for each model-scenario combination\n",
    "2. decadal/: Contains decade-aggregated files (2021-2030, 2031-2040, etc.)\n",
    "3. ensemble/: Contains multi-model ensemble statistics\n",
    "4. regional/: Contains country-level aggregated event counts\n",
    "5. ancillary/: Contains masks and threshold values\n",
    "\n",
    "File Naming Convention:\n",
    "----------------------\n",
    "GLEXIS_<model>_<scenario>_<year|decade>_<event>_<statistic>_v1.0.0.tif\n",
    "\n",
    "Variables:\n",
    "----------\n",
    "Each NetCDF file contains:\n",
    "- heat_days: Days with extreme heat\n",
    "- rain_days: Days with heavy precipitation\n",
    "- cold_days: Days with extreme cold\n",
    "- wind_days: Days with high wind speeds\n",
    "- drought_days: Legacy identifier for standardized water-balance lower-tail days\n",
    "- Various masks and threshold variables\n",
    "\n",
    "Models Included:\n",
    "---------------\n",
    paste("-", models, collapse = "\n"), "\n",
    "\n",
    "Scenarios Included:\n",
    "------------------\n",
    paste("-", scenarios, collapse = "\n"), "\n",
    "\n",
    "License:\n",
    "-------\n",
    "Creative Commons Attribution 4.0 International (CC BY 4.0)\n",
    "\n",
    "Citation:\n",
    "--------\n",
    pkg_metadata$reference, "\n",
    "DOI:", pkg_metadata$doi, "\n",
    "\n",
    "Contact:\n",
    "-------\n",
    pkg_metadata$author, "<", pkg_metadata$email, ">\n",
    pkg_metadata$institution, "\n",
    "\n",
    "Processing Date:\n",
    "---------------\n",
    format(Sys.time(), "%Y-%m-%d"), "\n"
  )
  
  writeLines(readme_content, file.path(output_dir, "documentation", "README.txt"))
  
  # Create data dictionary
  data_dict <- list(
    variables = list(
      heat_days = list(
        description = "Days with maximum temperature > 95th percentile",
        units = "days",
        threshold = "95th percentile of daily maximum temperature (1990-2010 baseline)"
      ),
      rain_days = list(
        description = "Heavy precipitation days > 95th percentile",
        units = "days",
        threshold = "95th percentile of daily precipitation (1990-2010 baseline)",
        note = ""
      ),
      cold_days = list(
        description = "Days with minimum temperature < 1st percentile",
        units = "days",
        threshold = "1st percentile of daily minimum temperature (1990-2010 baseline)"
      ),
      wind_days = list(
        description = "High wind speed days > 95th percentile",
        units = "days",
        threshold = "95th percentile of daily wind speed (1990-2010 baseline)",
        note = "For models with max wind data, 99th percentile is used"
      ),
      drought_days = list(
        description = "Legacy identifier for standardized water-balance lower-tail days; not canonical SPEI or drought severity",
        units = "days",
        threshold = "Annual within-cell z-score of precipitation minus PET below the 5th-percentile threshold-file layer"
      )
    ),
    masks = list(
      land_mask = "Land areas excluding Antarctica and glaciers (1=land, 0=ocean)",
      glacier_mask = "Glaciated areas (1=glacier, 0=non-glacier)",
      ocean_mask = "Ocean areas (1=ocean, 0=land)",
      equatorial_mask = "Equatorial region (-5° to 5° latitude)",
      tropical_mask = "Tropical region (-23.5° to 23.5° latitude)",
      arctic_mask = "Arctic region (>60°N latitude)",
      hyper_arid_mask = "Hyper-arid regions (January precipitation < 5 mm/month)",
      arid_mask = "Arid regions (annual precipitation < 30 mm/year)",
      permanent_wet_mask = "Permanently wet regions (annual precipitation > 2000 mm/year)",
      low_january_precip_mask = "Low January precipitation regions (<10 mm/month)"
    ),
    thresholds = list(
      tasmax_p95 = "95th percentile of daily maximum temperature (degrees C)",
      pr_p95 = "95th percentile of daily precipitation (mm/day)",
      tasmin_p01 = "1st percentile of daily minimum temperature (degrees C)",
      sfcwind_p95 = "95th percentile of daily wind speed (m/s)",
      water_balance_p05 = "5th-percentile threshold-file layer for the standardized water-balance proxy (unitless)"
    )
  )
  
  write_yaml(data_dict, file.path(output_dir, "documentation", "data_dictionary.yaml"))
  
  # Create processing manifest
  manifest <- list(
    processing = list(
      date = format(Sys.time(), "%Y-%m-%d"),
      parameters = params,
      models = models,
      scenarios = scenarios,
      years = paste(min(years), max(years), sep = "-"),
      targeted_fixes = params$targeted_fixes,
      software = list(
        name = pkg_metadata$name,
        version = pkg_metadata$version,
        r_version = R.version.string,
        packages = list(
          ncdf4 = as.character(packageVersion("ncdf4")),
          terra = as.character(packageVersion("terra")),
          SPEI = as.character(packageVersion("SPEI")),
          sf = as.character(packageVersion("sf")),
          rnaturalearth = as.character(packageVersion("rnaturalearth")),
          lubridate = as.character(packageVersion("lubridate")),
          yaml = as.character(packageVersion("yaml"))
        )
      )
    )
  )
  
  write_yaml(manifest, file.path(output_dir, "documentation", "processing_manifest.yaml"))
  
  cat("Documentation generated in:", file.path(output_dir, "documentation"), "\n")
}

# =============== MAIN PROCESSING FUNCTION ===============
process_single_year <- function(model, scenario, year, threshold_set = params$thresholds_quantiles, ensemble_index = NULL) {
  start_time <- Sys.time()
  cat("\n=== PROCESSING YEAR", year, "===\n")
  
  # Get model-specific pattern
  model_pattern <- get_model_pattern(model)
  
  # 1. FIND ALL REQUIRED FILES
  cat("Step 1: Finding all required files...\n")
  
  # Find which chunk contains our year
  chunk <- NULL
  for (chunk_name in names(file_chunks)) {
    if (year %in% file_chunks[[chunk_name]]) {
      chunk <- chunk_name
      break
    }
  }
  
  if (is.null(chunk)) {
    cat("Year", year, "not found in any file chunk\n")
    return(NULL)
  }
  
  # Get file paths for each variable
  file_paths <- list(
    tasmax = file.path(input_root, scenario, model,
                       paste0(model_pattern, "_", scenario, 
                              "_tasmax_global_daily_", chunk, ".nc")),
    pr = file.path(input_root, scenario, model,
                   paste0(model_pattern, "_", scenario, 
                          "_pr_global_daily_", chunk, ".nc")),
    tasmin = file.path(input_root, scenario, model,
                       paste0(model_pattern, "_", scenario, 
                              "_tasmin_global_daily_", chunk, ".nc")),
    sfcwind = file.path(input_root, scenario, model,
                        paste0(model_pattern, "_", scenario, 
                               "_sfcwind_global_daily_", chunk, ".nc")),
    hurs = file.path(input_root, scenario, model,
                     paste0(model_pattern, "_", scenario, 
                            "_hurs_global_daily_", chunk, ".nc")),
    ps = file.path(input_root, scenario, model,
                   paste0(model_pattern, "_", scenario, 
                          "_ps_global_daily_", chunk, ".nc")),
    pet = file.path(pet_root, scenario, model,
                    paste0(model_pattern, "_", scenario, 
                           "_pet_global_daily_", chunk, ".nc"))
  )
  
  # Check if all files exist
  missing_files <- sapply(file_paths, function(x) !file.exists(x))
  if (any(missing_files)) {
    cat("Missing files for year", year, ":\n")
    print(file_paths[missing_files])
    return(NULL)
  }
  
  # 2. GET GRID DIMENSIONS AND DATES
  cat("\nStep 2: Getting grid dimensions and dates...\n")
  
  # Use tasmax file as reference
  nc <- nc_open(file_paths$tasmax)
  lons <- ncvar_get(nc, "lon")
  lats <- ncvar_get(nc, "lat")
  time <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time")$units
  origin <- as.Date(gsub("days since ", "", time_units))
  dates <- as.Date(time, origin = origin)
  year_idx <- which(format(dates, "%Y") == sprintf("%04d", year))
  nc_close(nc)
  
  # Remove only South Pole (Antarctica)
  keep_lats <- which(lats >= params$south_pole_cutoff)
  lats <- lats[keep_lats]
  
  cat("Grid dimensions:", length(lons), "longitudes x", length(lats), "latitudes (Antarctica removed)\n")
  cat("Selected year:", year, ", Days:", length(year_idx), "\n")
  
  # 3. CREATE LAND MASK (excluding Antarctica and glaciers)
  cat("\nStep 3: Creating land mask...\n")
  world <- ne_countries(scale = 50, returnclass = "sf")
  glaciers <- ne_download(scale = 50, type = 'glaciated_areas', category = 'physical', returnclass = "sf")
  
  # Create base raster
  r <- rast(
    ncols = length(lons),
    nrows = length(lats),
    xmin = min(lons),
    xmax = max(lons),
    ymin = min(lats),
    ymax = max(lats),
    crs = "EPSG:4326"
  )
  
  # Rasterize land and exclude glaciers
  land_raster <- rasterize(world, r, field = 1, background = 0)
  glacier_raster <- rasterize(glaciers, r, field = 1, background = 0)
  land_mask <- land_raster * (1 - glacier_raster)
  land_mask_matrix <- t(as.matrix(land_mask, wide = TRUE))
  land_mask_3d <- array(land_mask_matrix, dim = c(dim(land_mask_matrix), length(year_idx)))
  land_mask_3d[land_mask_3d < 0.5] <- 0
  
  # 4. LOAD HISTORICAL CLIMATE DATA FOR DESERT/WET REGION DETECTION
  cat("\nStep 4: Loading historical climate data for desert/wet region detection...\n")
  
  # =============== PRECIPITATION THRESHOLD LOADING ===============
  load_pr_thresholds <- function(model) {
    # Define model-specific file patterns
    model_patterns <- list(
      "GFDL-ESM4" = "gfdl-esm4_r1i1p1f1_w5e5",
      "IPSL-CM6A-LR" = "ipsl-cm6a-lr_r1i1p1f1_w5e5",
      "MPI-ESM1-2-HR" = "mpi-esm1-2-hr_r1i1p1f1_w5e5",
      "MRI-ESM2-0" = "mri-esm2-0_r1i1p1f1_w5e5",
      "UKESM1-0-LL" = "ukesm1-0-ll_r1i1p1f2_w5e5"
    )
    
    # Get the correct pattern for this model
    pattern <- model_patterns[[model]]
    
    # Construct the threshold file path
    threshold_file <- file.path(thresholds_dir, 
                                paste0(pattern, "_thresholds_1990-2010_SPEI3.nc"))
    
    pr_thresholds <- ncvar_get(nc, "pr")[, , 4]  # Using 4th percentile (95th percentile)
    nc_close(nc)
    
    # Subset to exclude Antarctica
    pr_thresholds <- pr_thresholds[, keep_lats]
    
    # Calculate mean values (since threshold file doesn't contain time series)
    # These are approximations based on the thresholds
    mean_annual_precip <- pr_thresholds * 365  # Convert daily threshold to annual estimate
    mean_jan_precip <- pr_thresholds          # Using same values for January
    
    return(list(
      thresholds = pr_thresholds,
      mean_annual = mean_annual_precip,
      mean_jan = mean_jan_precip
    ))
  }
  
  # In your main processing, replace the precipitation loading with:
  precip_data <- tryCatch({
    pr_thresholds <- load_pr_thresholds(model)
    
    # Calculate derived values
    list(
      thresholds = pr_thresholds,
      mean_annual = pr_thresholds * 365,  # Convert daily to annual
      mean_jan = pr_thresholds            # Same for January
    )
  }, error = function(e) {
    # Fallback values (1 mm/day threshold)
    list(
      thresholds = matrix(1, nrow=length(lons), ncol=length(lats)),
      mean_annual = matrix(365, nrow=length(lons), ncol=length(lats)),
      mean_jan = matrix(30, nrow=length(lons), ncol=length(lats))
    )
  })
  
  
  if (!is.null(precip_data)) {
    mean_annual_precip <- precip_data$mean_annual
    mean_jan_precip <- precip_data$mean_jan
    cat("Mean annual precipitation range:", min(mean_annual_precip, na.rm = TRUE), 
        "-", max(mean_annual_precip, na.rm = TRUE), "mm/year\n")
    cat("Mean January precipitation range:", min(mean_jan_precip, na.rm = TRUE), 
        "-", max(mean_jan_precip, na.rm = TRUE), "mm/month\n")
  } else {
    mean_annual_precip <- NULL
    mean_jan_precip <- NULL
  }
  
  # Create hyper-arid mask (Jan precip < 5 mm/month)
  if (!is.null(mean_jan_precip)) {
    hyper_arid_mask <- mean_jan_precip < params$hyper_arid_threshold
    cat("Hyper-arid regions (Jan precip <", params$hyper_arid_threshold, "mm/month):", 
        sum(hyper_arid_mask, na.rm = TRUE), "pixels\n")
    hyper_arid_mask_3d <- array(hyper_arid_mask, dim = c(dim(hyper_arid_mask), length(year_idx)))
  } else {
    hyper_arid_mask <- matrix(FALSE, nrow = length(lons), ncol = length(lats))
    hyper_arid_mask_3d <- array(FALSE, dim = c(length(lons), length(lats), length(year_idx)))
  }
  
  # Create arid mask (annual precip < 30 mm/year)
  if (!is.null(mean_annual_precip)) {
    arid_mask <- mean_annual_precip < params$annual_arid_threshold
    cat("Arid regions (annual precip <", params$annual_arid_threshold, "mm/year):", 
        sum(arid_mask, na.rm = TRUE), "pixels\n")
    arid_mask_3d <- array(arid_mask, dim = c(dim(arid_mask), length(year_idx)))
  } else {
    arid_mask <- matrix(FALSE, nrow = length(lons), ncol = length(lats))
    arid_mask_3d <- array(FALSE, dim = c(length(lons), length(lats), length(year_idx)))
  }
  
  # Create permanent wet mask (annual precip > 2000 mm/year)
  if (!is.null(mean_annual_precip)) {
    permanent_wet_mask <- mean_annual_precip > params$annual_wet_threshold
    cat("Permanently wet regions (annual precip >", params$annual_wet_threshold, "mm/year):", 
        sum(permanent_wet_mask, na.rm = TRUE), "pixels\n")
    permanent_wet_mask_3d <- array(permanent_wet_mask, dim = c(dim(permanent_wet_mask), length(year_idx)))
  } else {
    permanent_wet_mask <- matrix(FALSE, nrow = length(lons), ncol = length(lats))
    permanent_wet_mask_3d <- array(FALSE, dim = c(length(lons), length(lats), length(year_idx)))
  }
  
  # Create January precipitation mask (Jan precip < 10 mm/month)
  if (!is.null(mean_jan_precip)) {
    jan_precip_mask <- mean_jan_precip < params$jan_precip_threshold
    cat("Low January precipitation regions (Jan precip <", params$jan_precip_threshold, "mm/month):", 
        sum(jan_precip_mask, na.rm = TRUE), "pixels\n")
    jan_precip_mask_3d <- array(jan_precip_mask, dim = c(dim(jan_precip_mask), length(year_idx)))
  } else {
    jan_precip_mask <- matrix(FALSE, nrow = length(lons), ncol = length(lats))
    jan_precip_mask_3d <- array(FALSE, dim = c(length(lons), length(lats), length(year_idx)))
  }
  
  # 5. CREATE HIGH-RESOLUTION OCEAN MASK
  cat("\nStep 5: Creating ocean mask...\n")
  coastline <- ne_coastline(scale = 50, returnclass = "sf")
  ocean_mask <- rasterize(coastline, r, field = 1, background = 0)
  ocean_mask_matrix <- t(as.matrix(ocean_mask, wide = TRUE))
  ocean_mask_3d <- array(ocean_mask_matrix, dim = dim(land_mask_3d))
  
  # 6. CREATE LATITUDE-BASED MASKS
  cat("\nStep 6: Creating latitude-based masks...\n")
  
  # Create latitude grid
  lat_grid <- array(rep(lats, each = length(lons)), dim = c(length(lons), length(lats)))
  
  # Equatorial ocean mask (-5° to 5° latitude)
  equatorial_ocean_mask <- (lat_grid >= -params$equatorial_lat_limit & lat_grid <= params$equatorial_lat_limit) & (ocean_mask_matrix == 1)
  equatorial_ocean_mask_3d <- array(equatorial_ocean_mask, dim = dim(land_mask_3d))
  
  # Tropical ocean mask (entire tropics)
  tropical_ocean_mask <- (lat_grid >= -params$tropical_lat_limit & lat_grid <= params$tropical_lat_limit) & (ocean_mask_matrix == 1)
  tropical_ocean_mask_3d <- array(tropical_ocean_mask, dim = dim(land_mask_3d))
  
  # Arctic mask (>60°N)
  arctic_mask <- lat_grid > 60
  arctic_mask_3d <- array(arctic_mask, dim = dim(land_mask_3d))
   
  # 7. CREATE COMPREHENSIVE EXCLUSION MASKS
  cat("\nStep 7: Creating comprehensive exclusion masks...\n")
  
  # Rain exclusion mask: Hyper-arid regions OR low January precipitation
  rain_exclude_mask <- hyper_arid_mask | jan_precip_mask
  cat("Rain exclusion regions:", sum(rain_exclude_mask, na.rm = TRUE), "pixels\n")
  
  # Cold exclusion mask: Ocean pixels
  cold_exclude_mask <- ocean_mask_matrix == 1
  cat("Cold exclusion regions:", sum(cold_exclude_mask, na.rm = TRUE), "pixels\n")
  
  # Heat exclusion mask: Equatorial ocean pixels
  heat_exclude_mask <- equatorial_ocean_mask
  cat("Heat exclusion regions:", sum(heat_exclude_mask, na.rm = TRUE), "pixels\n")
  
  # Wind exclusion mask: Ocean pixels
  wind_exclude_mask <- ocean_mask_matrix == 1
  cat("Wind exclusion regions:", sum(wind_exclude_mask, na.rm = TRUE), "pixels\n")
  
  # Wind Gust exclusion mask: Ocean pixels
  windgust_exclude_mask <- ocean_mask_matrix == 1
  cat("Wind Gust exclusion regions:", sum(windgust_exclude_mask, na.rm = TRUE), "pixels\n")
  
  # Drought exclusion mask: Arid regions OR permanently wet regions
  drought_exclude_mask <- arid_mask | permanent_wet_mask
  cat("Drought exclusion regions:", sum(drought_exclude_mask, na.rm = TRUE), "pixels\n")
  
  # Create 3D versions of exclusion masks
  rain_exclude_mask_3d <- array(rain_exclude_mask, dim = c(dim(rain_exclude_mask), length(year_idx)))
  cold_exclude_mask_3d <- array(cold_exclude_mask, dim = c(dim(cold_exclude_mask), length(year_idx)))
  heat_exclude_mask_3d <- array(heat_exclude_mask, dim = c(dim(heat_exclude_mask), length(year_idx)))
  wind_exclude_mask_3d <- array(wind_exclude_mask, dim = c(dim(wind_exclude_mask), length(year_idx)))
  windgust_exclude_mask_3d <- array(windgust_exclude_mask, dim = c(dim(windgust_exclude_mask), length(year_idx)))
  drought_exclude_mask_3d <- array(drought_exclude_mask, dim = c(dim(drought_exclude_mask), length(year_idx)))
  

  # 7. LOAD THRESHOLDS - FINAL WORKING VERSION
  cat("\nStep 7: Loading thresholds...\n")
  
  # Define threshold quantiles (original plus windgust)
  threshold_quantiles <- list(
    tasmax = 4,      # 95th percentile
    pr = 4,          # 95th percentile
    tasmin = 1,      # p01 threshold-file layer
    sfcwind = 4,     # 95th percentile
    windgust = 4,    # 95th percentile
    spei = 2         # p05 baseline SPEI-3 threshold-file layer
  )
  
  # Function to load thresholds with exact file matching
  load_thresholds <- function(var) {
    # Define model base name without run identifier
    model_base <- tolower(gsub("-", "_", model))
    
    # SPECIAL HANDLING FOR WIND VARIABLES
    if (var %in% c("sfcwind", "windgust")) {
      # Use ONLY wind-specific directory
      possible_files <- c(
        file.path(wind_thresholds_dir, paste0(model_base, "_", var, "_thresholds_1990-2010.nc"))
      )
      
      # Find first existing file
      threshold_file <- NULL
      for (f in possible_files) {
        if (file.exists(f)) {
          threshold_file <- f
          break
        }
      }
      
      if (is.null(threshold_file)) {
        stop(paste("No threshold file found for", var, 
                   "\nSearched for:\n", paste(possible_files, collapse = "\n")))
      }
      
      cat("Loading thresholds for", var, "from:", threshold_file, "\n")
      nc <- nc_open(threshold_file)
      
      # Handle variable names
      var_names <- names(nc$var)
      target_var <- ifelse(var == "sfcwind", "sfcWind", "windGust")
      
      if (!target_var %in% var_names) {
        nc_close(nc)
        stop(paste("Variable", target_var, "not found in", threshold_file))
      }
      
      data <- ncvar_get(nc, target_var)
      nc_close(nc)
      
      # Subset to exclude Antarctica
      data <- data[, keep_lats, ]
      
      return(data)
    }
    
    # ORIGINAL CODE FOR NON-WIND VARIABLES
    # Use the original model pattern with dashes
    model_pattern <- tolower(model)
    
    # Define base pattern with run identifier
    base_pattern <- paste0(model_pattern, "_r1i1p1f1_w5e5")
    
    # Define possible file locations and patterns
    possible_files <- c(
      file.path(thresholds_dir, paste0(base_pattern, "_thresholds_1990-2010_SPEI3.nc")),
      file.path(thresholds_dir, paste0(base_pattern, "_", var, "_thresholds_1990-2010.nc")),
      file.path(thresholds_dir, paste0(base_pattern, "_thresholds_1990-2010_", var, ".nc")),
      # Fallback patterns without the run identifier
      file.path(thresholds_dir, paste0(model_pattern, "_thresholds_1990-2010_SPEI3.nc")),
      file.path(thresholds_dir, paste0(model_pattern, "_", var, "_thresholds_1990-2010.nc")),
      file.path(thresholds_dir, paste0(model_pattern, "_thresholds_1990-2010_", var, ".nc"))
    )
    
    # Find first existing file
    threshold_file <- NULL
    for (f in possible_files) {
      if (file.exists(f)) {
        threshold_file <- f
        break
      }
    }
    
    if (is.null(threshold_file)) {
      stop(paste("No threshold file found for", var, 
                 "\nSearched for:\n", paste(possible_files, collapse = "\n")))
    }
    
    cat("Loading thresholds for", var, "from:", threshold_file, "\n")
    nc <- nc_open(threshold_file)
    
    # Handle different variable names in files
    var_names <- names(nc$var)
    target_var <- var
    if (!target_var %in% var_names) {
      target_var <- var_names[grep(var, var_names)[1]]
    }
    
    data <- ncvar_get(nc, target_var)
    nc_close(nc)
    
    # Subset to exclude Antarctica
    data <- data[, keep_lats, ]
    
    return(data)
  }
  
  # Load thresholds for each variable
  thresholds <- list()
  for (var in c("tasmax", "pr", "tasmin", "sfcwind", "windgust", "spei")) {
    tryCatch({
      thresholds[[var]] <- load_thresholds(var)
    }, error = function(e) {
      cat("Error loading thresholds for", var, ":", e$message, "\n")
      stop("Threshold loading failed")
    })
  }
  # 8. CREATE 3D ANNUAL THRESHOLD ARRAYS
  cat("\nStep 8: Creating 3D threshold arrays...\n")
  thresholds_3d <- list()
  
  for (var in names(thresholds)) {
    thresh_data <- thresholds[[var]]
    dims <- dim(thresh_data)
    
    cat("Processing", var, "with dimensions:", paste(dims, collapse = "x"), "\n")
    
    # Special handling for SPEI thresholds
    if (var == "spei") {
      # SPEI thresholds are typically constant values
      if (length(dims) == 3 && dims[3] == 9) {
        # Select the appropriate quantile (2nd for SPEI)
        quantile_slice <- thresh_data[, , threshold_quantiles[[var]]]
        
        # Create 3D array by replicating for each day
        thresholds_3d[[var]] <- array(quantile_slice, 
                                      dim = c(dim(quantile_slice), length(year_idx)))
      } else {
        stop(paste("Unexpected SPEI threshold dimensions:", paste(dims, collapse = "x")))
      }
    }
    # Case 1: If we have multiple quantiles stored in different dimensions
    else if (length(dims) == 4) {
      # Select the appropriate quantile
      thresholds_3d[[var]] <- thresh_data[, , , threshold_quantiles[[var]]]
    } 
    # Case 2: If we have a 3D array with 5 layers (likely quantiles)
    else if (length(dims) == 3 && dims[3] == 5) {
      # Select the appropriate quantile (assuming 3rd dimension is quantiles)
      quantile_slice <- thresh_data[, , threshold_quantiles[[var]]]
      
      # Create 3D array by replicating for each day
      thresholds_3d[[var]] <- array(quantile_slice, 
                                    dim = c(dim(quantile_slice), length(year_idx)))
    }
    # Case 3: Daily climatology (365 days)
    else if (length(dims) == 3 && dims[3] == 365) {
      # Select days matching our year
      thresholds_3d[[var]] <- thresh_data[, , yday(dates[year_idx])]
    }
    # Case 4: Monthly thresholds (12 months)
    else if (length(dims) == 3 && dims[3] == 12) {
      # Replicate monthly thresholds for each day in our year
      thresholds_3d[[var]] <- array(NA, dim = c(dims[1:2], length(year_idx)))
      for (m in 1:12) {
        month_days <- which(month(dates[year_idx]) == m)
        if (length(month_days) > 0) {
          thresholds_3d[[var]][, , month_days] <- array(thresh_data[, , m], 
                                                        dim = c(dims[1:2], length(month_days)))
        }
      }
    }
    # Case 5: Single 2D threshold
    else if (length(dims) == 2) {
      # Replicate the threshold for each day
      thresholds_3d[[var]] <- array(thresh_data, 
                                    dim = c(dims, length(year_idx)))
    }
    else {
      stop(paste("Unexpected threshold dimensions for", var, 
                 ":", paste(dims, collapse = "x")))
    }
    
    # Final dimension check
    if (length(dim(thresholds_3d[[var]])) != 3) {
      stop(paste("Failed to create proper 3D array for", var,
                 "Resulting dimensions:", 
                 paste(dim(thresholds_3d[[var]]), collapse = "x")))
    }
    
    # Ensure land mask has the same dimensions
    if (!identical(dim(thresholds_3d[[var]]), dim(land_mask_3d))) {
      cat("Reshaping land mask to match threshold dimensions\n")
      land_mask_reshaped <- array(land_mask_3d[, , 1], dim = dim(thresholds_3d[[var]]))
    } else {
      land_mask_reshaped <- land_mask_3d
    }
    
    # Apply land mask
    thresholds_3d[[var]] <- thresholds_3d[[var]] * land_mask_reshaped
  }
  
  # 9. LOAD CLIMATE DATA (with wind max preference)
  cat("\nStep 9: Loading climate data...\n")
  
  # Track wind data type
  wind_data_type <- "mean"  # Default to mean wind
  
  load_var <- function(var, unit_conversion = NULL) {
    # Special handling for wind: prefer max data if available
    if (var == "sfcwind" || var == "windgust") {
      # First try with run identifier
      max_file1 <- file.path(input_root, scenario, model,
                             paste0(tolower(model), "_r1i1p1f1_w5e5_", scenario, 
                                    "_", var, "max_global_daily_", chunk, ".nc"))
      max_file2 <- file.path(input_root, scenario, model,
                             paste0(tolower(model), "_", scenario, 
                                    "_", var, "max_global_daily_", chunk, ".nc"))
      max_file3 <- file.path(input_root, scenario, model,
                             paste0(model_pattern, "_", scenario, 
                                    "_", var, "_global_daily_", chunk, ".nc"))
      max_file4 <- file.path(input_root, scenario, model,
                             paste0("ukesm1-0-ll_r1i1p1f2_w5e5_", scenario, 
                                    "_", var, "_global_daily_", chunk, ".nc"))
      
      if (file.exists(max_file1)) {
        cat("Using", var, "MAX data (with run identifier)\n")
        nc_file <- max_file1
        var_name <- paste0(var, "max")
        wind_data_type <<- "max"
      } else if (file.exists(max_file2)) {
        cat("Using", var, "MAX data (without run identifier)\n")
        nc_file <- max_file2
        var_name <- paste0(var, "max")
        wind_data_type <<- "max"
      } else {
        # Try mean wind data
        mean_file1 <- file.path(input_root, scenario, model,
                                paste0(tolower(model), "_r1i1p1f1_w5e5_", scenario, 
                                       "_", var, "_global_daily_", chunk, ".nc"))
        mean_file2 <- file.path(input_root, scenario, model,
                                paste0(tolower(model), "_", scenario, 
                                       "_", var, "_global_daily_", chunk, ".nc"))
        mean_file3 <- file.path(input_root, scenario, model,
                                paste0(model_pattern, "_", scenario, 
                                       "_", var, "_global_daily_", chunk, ".nc"))
        mean_file4 <- file.path(input_root, scenario, model,
                                paste0("ukesm1-0-ll_r1i1p1f2_w5e5_", scenario, 
                                       "_", var, "_global_daily_", chunk, ".nc"))
        
        if (file.exists(mean_file1)) {
          cat("Using", var, "MEAN data (with run identifier)\n")
          nc_file <- mean_file1
          var_name <- var
          wind_data_type <<- "mean"
        } else if (file.exists(mean_file2)) {
          cat("Using", var, "MEAN data (without run identifier)\n")
          nc_file <- mean_file2
          var_name <- var
          wind_data_type <<- "mean"
        } else if (file.exists(mean_file3)) {
          cat("Using", var, "MEAN data (ukesm)\n")
          nc_file <- mean_file3
          var_name <- var
          wind_data_type <<- "mean"
        } else if (file.exists(mean_file4)) {
          cat("Using", var, "MEAN data (ukesm)\n")
          nc_file <- mean_file4
          var_name <- var
          wind_data_type <<- "mean"
        } else {
          stop(paste(var, "data file not found. Tried:\n", 
                     max_file1, "\n", max_file2, "\n", 
                     mean_file1, "\n", mean_file2))
        }
      }
    } else {
      # For non-wind variables
      file1 <- file.path(input_root, scenario, model,
                         paste0(tolower(model), "_r1i1p1f1_w5e5_", scenario, 
                                "_", var, "_global_daily_", chunk, ".nc"))
      file2 <- file.path(input_root, scenario, model,
                         paste0(tolower(model), "_", scenario, 
                                "_", var, "_global_daily_", chunk, ".nc"))
      file3 <- file.path(input_root, scenario, model,
                         paste0(model_pattern, "_", scenario, 
                                "_", var, "_global_daily_", chunk, ".nc"))
      file4 <- file.path(input_root, scenario, model,
                         paste0("ukesm1-0-ll_r1i1p1f2_w5e5_", scenario, 
                                "_", var, "_global_daily_", chunk, ".nc"))
      
      if (file.exists(file1)) {
        nc_file <- file1
        cat("Using file with run identifier\n")
      } else if (file.exists(file2)) {
        nc_file <- file2
        cat("Using file without run identifier\n")
      } else if (file.exists(file3)) {
        nc_file <- file3
        cat("Using file for ukesm\n")
      } else if (file.exists(file4)) {
        nc_file <- file4
        cat("Using file for ukesm\n")
      } else {
        stop(paste("Climate data file not found for", var, 
                   "\nTried:\n", file1, "\n", file2))
      }
      var_name <- var
    }
    
    if (!file.exists(nc_file)) {
      stop(paste("Climate data file not found:", nc_file))
    }
    
    cat("Loading:", nc_file, "\n")
    nc <- nc_open(nc_file)
    
    # Read full grid then subset
    data <- ncvar_get(nc, var_name, 
                      start = c(1, 1, min(year_idx)), 
                      count = c(-1, -1, length(year_idx)))
    nc_close(nc)
    
    # Subset to exclude Antarctica
    data <- data[, keep_lats, ]
    
    if (!is.null(unit_conversion)) {
      if (unit_conversion == "K_to_C") data <- data - 273.15
      if (unit_conversion == "kgm2s_to_mmday") data <- data * 86400
    }
    
    return(data)
  }  
  
  # Load data with units conversion
  year_data <- list(
    tasmax = load_var("tasmax", "K_to_C"),
    pr = load_var("pr", "kgm2s_to_mmday"),
    tasmin = load_var("tasmin", "K_to_C"),
    sfcwind = load_var("sfcwind"),
    hurs = load_var("hurs"),
    ps = load_var("ps")
  )
  
  # If windgust data not available, estimate it from other variables
  if (!exists("windgust", where = year_data)) {
    cat("Estimating wind gust from available data...\n")
    year_data$windgust <- estimate_wind_gust_vectorized(
      mean_wind = year_data$sfcwind,
      temperature = year_data$tasmin,
      humidity = year_data$hurs,
      pressure = year_data$ps
    )
  }
  
  # Load PET data
  pet_file1 <- file.path(pet_root, scenario, model,
                         paste0(tolower(model), "_r1i1p1f1_w5e5_", scenario,
                                "_pet_global_daily_", chunk, ".nc"))
  pet_file2 <- file.path(pet_root, scenario, model,
                         paste0(tolower(model), "_", scenario,
                                "_pet_global_daily_", chunk, ".nc"))
  pet_file3 <- file.path(pet_root, scenario, model,
                         paste0(tolower(model), "_r1i1p1f2_w5e5_", scenario,
                                "_pet_global_daily_", chunk, ".nc"))
  
  if (file.exists(pet_file1)) {
    pet_file <- pet_file1
  } else if (file.exists(pet_file2)) {
    pet_file <- pet_file2
  } else {
    pet_file <- pet_file3
  }
  if (!file.exists(pet_file)) {
    stop(paste("PET file not found. Tried:\n", pet_file1, "\n", pet_file2))
  }
  
  cat("Loading PET:", pet_file, "\n")
  nc_pet <- nc_open(pet_file)

  pet_units <- ncatt_get(nc_pet, "pet", "units")$value
  if (is.null(pet_units) || !tolower(trimws(pet_units)) %in% c("mm/day", "mm d-1", "mm day-1")) {
    nc_close(nc_pet)
    stop("PET must be supplied in daily water-equivalent units (mm/day); found: ", pet_units)
  }
  
  # Get time information
  time_pet <- ncvar_get(nc_pet, "time")
  time_units <- ncatt_get(nc_pet, "time")$units
  calendar <- ncatt_get(nc_pet, "time", "calendar")$value
  if (is.null(calendar)) calendar <- "standard"
  
  origin <- as.Date(gsub("days since ", "", time_units))
  dates_pet <- as.Date(time_pet, origin = origin)
  
  # Find available years in the file
  file_years <- unique(format(dates_pet, "%Y"))
  
  # Verify requested year exists in file
  if (!as.character(year) %in% file_years) {
    nc_close(nc_pet)
    stop(paste("Year", year, "not found in PET file. Available years:", 
               paste(file_years, collapse = ", ")))
  }
  
  # Find indices for our target year
  year_idx <- which(format(dates_pet, "%Y") == sprintf("%04d", year))
  
  # Verify we found the year
  if (length(year_idx) == 0) {
    nc_close(nc_pet)
    stop(paste("No time indices found for year", year, "in PET file"))
  }
  
  # Read the data
  year_data$pet <- tryCatch({
    ncvar_get(nc_pet, "pet", 
              start = c(1, 1, min(year_idx)), 
              count = c(-1, -1, length(year_idx)))
  }, error = function(e) {
    nc_close(nc_pet)
    stop(paste("Failed to read PET data:", e$message))
  })
  
  nc_close(nc_pet)
  
  # Verify dimensions
  if (length(dim(year_data$pet)) != 3) {
    stop("PET data is not 3-dimensional")
  }
  
  # Subset to exclude Antarctica
  year_data$pet <- year_data$pet[, keep_lats, ]
  
  # 10. APPLY LAND MASK TO DATA
  cat("\nStep 10: Applying land mask...\n")
  
  time_dims <- sapply(year_data, function(x) dim(x)[3])
  ref_time_dim <- as.numeric(names(which.max(table(time_dims))))
  
  for (var in names(year_data)) {
    current_dims <- dim(year_data[[var]])
    
    if (current_dims[3] != ref_time_dim) {
      if (current_dims[3] > ref_time_dim) {
        year_data[[var]] <- year_data[[var]][, , -60]
      }
      else {
        fb <- year_data[[var]][, , 59]
        year_data[[var]] <- abind::abind(
          year_data[[var]][, , 1:59],
          fb,
          year_data[[var]][, , 60:365],
          along = 3
        )
      }
    }
  }
  
  if (dim(land_mask_3d)[3] != ref_time_dim) {
    if (dim(land_mask_3d)[3] > ref_time_dim) {
      land_mask_3d <- land_mask_3d[, , -60]  
    } else {
      feb28_mask <- land_mask_3d[, , 59]
      land_mask_3d <- abind::abind(
        land_mask_3d[, , 1:59],
        feb28_mask,
        land_mask_3d[, , 60:365],
        along = 3
      )
    }
  }
  
  for (var in names(year_data)) {
    if (!identical(dim(year_data[[var]]), dim(land_mask_3d))) {
      stop(paste("Dimension mismatch after alignment for", var,
                 "Expected:", paste(dim(land_mask_3d), collapse="x"),
                 "Got:", paste(dim(year_data[[var]]), collapse="x")))
    }
  }
  
  for (var in names(year_data)) {
    year_data[[var]] <- year_data[[var]] * land_mask_3d
  }
  
  # 12. CALCULATE WATER BALANCE AND ANNUAL WITHIN-CELL STANDARDIZED SCORES
  cat("\nStep 12: Calculating water balance and annual within-cell standardized scores...\n")
  
  # Calculate water balance (P - PET)
  wb <- year_data$pr - year_data$pet
  
  # Create SPEI result array
  spei_result <- array(NA, dim = dim(wb))
  
  # Get land points
  land_points <- which(land_mask_matrix == 1, arr.ind = TRUE)
  
  # Calculate SPEI using standardized anomalies
  for (i in 1:nrow(land_points)) {
    idx <- land_points[i, ]
    
    # Skip if outside data bounds
    if (idx[1] > dim(wb)[1] || idx[2] > dim(wb)[2]) next
    
    wb_series <- wb[idx[1], idx[2], ]
    
    if (sum(!is.na(wb_series)) >= params$spei_scale) {
      tryCatch({
        # Standardize using monthly climatology
        wb_mean <- mean(wb_series, na.rm = TRUE)
        wb_sd <- sd(wb_series, na.rm = TRUE)
        
        # Zero or non-finite annual standard deviations remain missing.
        if (is.finite(wb_sd) && wb_sd > 0) {
          spei_result[idx[1], idx[2], ] <- (wb_series - wb_mean) / wb_sd
        }
      }, error = function(e) {
        cat("SPEI error at", idx[1], idx[2], ":", e$message, "\n")
      })
    }
  }
  
  # Apply land mask and clean up
  spei_result <- spei_result * land_mask_3d
  spei_result[is.infinite(spei_result)] <- NA
  spei_result <- round(spei_result, 2)
  
  # 13. CALCULATE EXTREME EVENTS WITH COMPREHENSIVE EXCLUSIONS
  cat("\nStep 13: Calculating extreme events with comprehensive exclusions...\n")
  
  events <- list(
    # Heat: Use 95th percentile threshold but exclude equatorial oceans
    heat = (year_data$tasmax > thresholds_3d$tasmax) & 
      !equatorial_ocean_mask_3d,
    
    # Rain: Strict masking with minimum precipitation
    rain = (year_data$pr > thresholds_3d$pr) & 
      (year_data$pr >= params$min_precip_filter) &
      !hyper_arid_mask_3d,
    
    # Cold: Strict land mask only
    cold = (year_data$tasmin < thresholds_3d$tasmin) & 
      !equatorial_ocean_mask_3d,
    
    # Wind: Strict land mask and minimum speed for mean wind
    wind = (year_data$sfcwind > thresholds_3d$sfcwind) & 
      (land_mask_3d == 1),
    
    windgust = (year_data$windgust > thresholds_3d$windgust) & 
      (land_mask_3d == 1),
    
    # Drought: Mask permanently dry and wet regions
    drought = (spei_result < thresholds_3d$spei) & 
      !arid_mask_3d & 
      !permanent_wet_mask_3d
  )
  
  # Additional wind fix for mean wind data
  if (wind_data_type == "mean") {
    events$wind <- events$wind & (year_data$sfcwind >= params$wind_min_speed)
    events$windgust <- events$windgust & (year_data$windgust >= params$wind_min_speed)
  }
  
  # FINAL TARGETED FIXES FOR PROBLEMATIC PIXELS
  cat("\nApplying targeted fixes for remaining problematic pixels...\n")
  
  # 1. RAINFALL: Atacama Desert cells in January
  if (1 %in% month(dates[year_idx])) {
    jan_days <- which(month(dates[year_idx]) == 1)
    
    # Coordinates of problematic Atacama cells
    atacama_cells <- list(
      c(lon = -74.64618, lat = -14.65083),
      c(lon = -68.65451, lat = -28.60417)
    )
    
    for (cell in atacama_cells) {
      i <- which.min(abs(lons - cell["lon"]))
      j <- which.min(abs(lats - cell["lat"]))
      
      if (i <= dim(events$rain)[1] && j <= dim(events$rain)[2]) {
        events$rain[i, j, jan_days] <- FALSE
        cat("Zeroed out Atacama rainfall events at", cell["lon"], cell["lat"], "in January\n")
      }
    }
  }
  
  # 2. COLD: Caribbean ocean cell in February
  if (2 %in% month(dates[year_idx])) {
    feb_days <- which(month(dates[year_idx]) == 2)
    
    caribbean_cell <- c(lon = -66.65729, lat = 18.23917)
    i <- which.min(abs(lons - caribbean_cell["lon"]))
    j <- which.min(abs(lats - caribbean_cell["lat"]))
    
    if (i <= dim(events$cold)[1] && j <= dim(events$cold)[2]) {
      events$cold[i, j, feb_days] <- FALSE
      cat("Zeroed out Caribbean cold events at", caribbean_cell["lon"], caribbean_cell["lat"], "in February\n")
    }
  }
  
  # 3. WIND: North Atlantic ocean cell in September
  if (9 %in% month(dates[year_idx])) {
    sep_days <- which(month(dates[year_idx]) == 9)
    
    atlantic_cell <- c(lon = -70.65174, lat = 66.5775)
    i <- which.min(abs(lons - atlantic_cell["lon"]))
    j <- which.min(abs(lats - atlantic_cell["lat"]))
    
    if (i <= dim(events$wind)[1] && j <= dim(events$wind)[2]) {
      events$wind[i, j, sep_days] <- FALSE
      cat("Zeroed out North Atlantic wind events at", atlantic_cell["lon"], atlantic_cell["lat"], "in September\n")
    }
  }
  
  # 4. WINDGUST: North Atlantic ocean cell in September
  if (9 %in% month(dates[year_idx])) {
    sep_days <- which(month(dates[year_idx]) == 9)
    
    atlantic_cell <- c(lon = -70.65174, lat = 66.5775)
    i <- which.min(abs(lons - atlantic_cell["lon"]))
    j <- which.min(abs(lats - atlantic_cell["lat"]))
    
    if (i <= dim(events$windgust)[1] && j <= dim(events$windgust)[2]) {
      events$windgust[i, j, sep_days] <- FALSE
      cat("Zeroed out North Atlantic wind gust events at", atlantic_cell["lon"], atlantic_cell["lat"], "in September\n")
    }
  }
  
  # 5. WATER-BALANCE PROXY: targeted correction cell in June
  if (6 %in% month(dates[year_idx])) {
    jun_days <- which(month(dates[year_idx]) == 6)
    
    targeted_drought_cell <- c(lon = 50.67951, lat = 10.76417)
    i <- which.min(abs(lons - targeted_drought_cell["lon"]))
    j <- which.min(abs(lats - targeted_drought_cell["lat"]))
    
    if (i <= dim(events$drought)[1] && j <= dim(events$drought)[2]) {
      events$drought[i, j, jun_days] <- FALSE
      cat("Applied targeted drought correction at", targeted_drought_cell["lon"], targeted_drought_cell["lat"], "in June\n")
    }
  }
  
  # FINAL LAND MASK APPLICATION
  for (event in names(events)) {
    events[[event]] <- ifelse(land_mask_3d == 1, events[[event]], NA)
  }
  
  # 14. SUM EVENTS BY MONTH AND CALCULATE MAGNITUDE STATISTICS
  cat("\nStep 14: Summing events and calculating magnitude statistics...\n")
  
  # Initialize monthly event counts and magnitude arrays
  monthly_event_counts <- list()
  monthly_magnitude_stats <- list()
  
  for (event in names(events)) {
    monthly_event_counts[[event]] <- array(0, dim = c(dim(events[[event]])[1:2], 12))
    
    # Create structure for magnitude statistics (max, min, mean)
    monthly_magnitude_stats[[event]] <- array(NA, 
                                              dim = c(dim(events[[event]])[1:2], 12, 3),
                                              dimnames = list(NULL, NULL, NULL, c("max", "min", "mean"))
    )
  }
  
  # Calculate monthly sums and magnitude statistics
  for (m in 1:12) {
    month_days <- which(month(dates[year_idx]) == m)
    if (length(month_days) > 0) {
      for (event in names(events)) {
        # Count events
        event_slice <- events[[event]][, , month_days, drop = FALSE]
        valid_cells <- apply(!is.na(event_slice), c(1,2), any)
        month_counts <- apply(event_slice, c(1,2), sum, na.rm = TRUE)
        month_counts[!valid_cells] <- NA_real_
        monthly_event_counts[[event]][, , m] <- month_counts
        
        # Get corresponding data for magnitude calculations
        event_data <- switch(event,
                             "heat" = year_data$tasmax,
                             "rain" = year_data$pr,
                             "cold" = year_data$tasmin,
                             "wind" = year_data$sfcwind,      
                             "windgust" = year_data$windgust,  
                             "drought" = spei_result)
        
        # Apply event mask to get only extreme values
        masked_data <- event_data[, , month_days] * events[[event]][, , month_days]
        
        # Calculate magnitude statistics
        monthly_magnitude_stats[[event]][, , m, "max"] <- apply(masked_data, c(1,2), max, na.rm = TRUE)
        monthly_magnitude_stats[[event]][, , m, "min"] <- apply(masked_data, c(1,2), min, na.rm = TRUE)
        monthly_magnitude_stats[[event]][, , m, "mean"] <- apply(masked_data, c(1,2), mean, na.rm = TRUE)
        
        # Clean up infinite values
        monthly_magnitude_stats[[event]][, , m, ][is.infinite(monthly_magnitude_stats[[event]][, , m, ])] <- NA
      }
    }
  }
  
  # 15. CREATE OUTPUT RASTERS FOR COUNTS AND MAGNITUDE STATISTICS
  cat("\nStep 15: Creating output rasters...\n")
  
  create_output_raster <- function(data) {
    r <- rast(t(data), crs = "EPSG:4326")
    dx <- median(diff(sort(unique(lons))))
    dy <- median(diff(sort(unique(lats))))
    ext(r) <- c(min(lons) - dx / 2, max(lons) + dx / 2,
                min(lats) - dy / 2, max(lats) + dy / 2)
    return(r)
  }
  
  # Create monthly rasters for each event type (counts and magnitude stats)
  monthly_count_rasters <- list()
  monthly_magnitude_rasters <- list()
  
  for (event in names(monthly_event_counts)) {
    monthly_count_rasters[[event]] <- lapply(1:12, function(m) {
      create_output_raster(monthly_event_counts[[event]][, , m])
    })
    
    # Create separate rasters for each magnitude statistic
    monthly_magnitude_rasters[[event]] <- list(
      max = lapply(1:12, function(m) create_output_raster(monthly_magnitude_stats[[event]][, , m, "max"])),
      min = lapply(1:12, function(m) create_output_raster(monthly_magnitude_stats[[event]][, , m, "min"])),
      mean = lapply(1:12, function(m) create_output_raster(monthly_magnitude_stats[[event]][, , m, "mean"]))
    )
  }
  
  # 16. SAVE OUTPUTS WITH MAGNITUDE STATISTICS
  cat("\nStep 16: Saving outputs with magnitude statistics...\n")
  
  # Create output directory for this year
  year_dir <- file.path(output_dir, "annual", paste0(tolower(model), "_", scenario), year)
  dir.create(year_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create ensemble suffix if needed
  suffix <- ifelse(!is.null(ensemble_index), paste0("_ensemble", ensemble_index), "")
  
  # Save each event type with count and magnitude statistics
  for (event in names(monthly_count_rasters)) {
    # Save count data
    r_stack_count <- rast(monthly_count_rasters[[event]])
    names(r_stack_count) <- month.name[1:12]
    time(r_stack_count) <- as.Date(paste(year, 1:12, "15", sep = "-"))
    varnames(r_stack_count) <- "event_days"
    units(r_stack_count) <- "days"
    
    out_file_count <- file.path(year_dir, 
                                paste0("GLEXIS_", tolower(model), "_", scenario, "_", year,
                                       "_", event, suffix, "_monthly_counts.tif"))
    writeRaster(r_stack_count, filename = out_file_count, overwrite = TRUE,
                gdal = c("DESCRIPTION=Monthly extreme event counts"))
    cat("Saved:", out_file_count, "\n")
    
    # Save magnitude statistics (max, min, mean)
    for (stat in c("max", "min", "mean")) {
      r_stack_mag <- rast(monthly_magnitude_rasters[[event]][[stat]])
      names(r_stack_mag) <- month.name[1:12]
      time(r_stack_mag) <- as.Date(paste(year, 1:12, "15", sep = "-"))
      
      # Set variable-specific metadata
      if (event == "heat") {
        varnames(r_stack_mag) <- paste0(stat, "_tasmax")
        units(r_stack_mag) <- "°C"
        desc <- paste("Monthly", stat, "of daily maximum temperature during heat events")
      } else if (event == "rain") {
        varnames(r_stack_mag) <- paste0(stat, "_pr")
        units(r_stack_mag) <- "mm/day"
        desc <- paste("Monthly", stat, "of daily precipitation during heavy rain events")
      } else if (event == "cold") {
        varnames(r_stack_mag) <- paste0(stat, "_tasmin")
        units(r_stack_mag) <- "°C"
        desc <- paste("Monthly", stat, "of daily minimum temperature during cold spells")
      } else if (event == "wind") {
        varnames(r_stack_mag) <- paste0(stat, "_sfcwind")
        units(r_stack_mag) <- "m/s"
        desc <- paste("Monthly", stat, "of daily surface wind speed during wind events")
      } else if (event == "windgust") {
        varnames(r_stack_mag) <- paste0(stat, "_windgust")
        units(r_stack_mag) <- "m/s"
        desc <- paste("Monthly", stat, "of estimated daily wind gust during wind events")
      } else if (event == "drought") {
        varnames(r_stack_mag) <- paste0(stat, "_spei")
        units(r_stack_mag) <- "index"
        desc <- paste("Monthly", stat, "of SPEI values during drought events")
      }
      
      out_file_mag <- file.path(year_dir, 
                                paste0("GLEXIS_", tolower(model), "_", scenario, "_", year,
                                       "_", event, suffix, "_monthly_", stat, "_magnitude.tif"))
      writeRaster(r_stack_mag, filename = out_file_mag, overwrite = TRUE,
                  gdal = c(paste0("DESCRIPTION=", desc)))
      cat("Saved:", out_file_mag, "\n")
    }
  }
}
# =============== MAIN EXECUTION ===============
cat("GLEXIS Processing Pipeline\n")
cat("=========================\n")
cat("Version:", pkg_metadata$version, "\n")
cat("Author:", pkg_metadata$author, "\n")
cat("Institution:", pkg_metadata$institution, "\n")
cat("License:", pkg_metadata$license, "\n")
cat("DOI:", pkg_metadata$doi, "\n\n")

global_start_time <- Sys.time()

# Initialize storage for ensemble processing
ensemble_data <- list()

# Process based on run type
if (params$run_type == "main") {
  # Process each year with main thresholds
  for (year in years) {
    result <- tryCatch({
      process_single_year(model, scenario, year)
    }, error = function(e) {
      cat("Error processing", model, scenario, year, ":", e$message, "\n")
      NULL
    })
  }
} else if (params$run_type == "ensemble") {
  # Process with each ensemble threshold set
  n_ensemble <- length(params$ensemble_thresholds$tasmax)
  for (i in 1:n_ensemble) {
    cat("\n===== PROCESSING ENSEMBLE SET", i, "=====\n")
    ensemble_threshold_set <- list(
      tasmax = params$ensemble_thresholds$tasmax[i],
      pr = params$ensemble_thresholds$pr[i],
      tasmin = params$ensemble_thresholds$tasmin[i],
      sfcwind = params$ensemble_thresholds$sfcwind[i],
      windgust = params$ensemble_thresholds$windgust[i],
      spei = params$ensemble_thresholds$spei[i]
    )
    
    for (year in years) {
      result <- tryCatch({
        process_single_year(model, scenario, year, 
                            threshold_set = ensemble_threshold_set,
                            ensemble_index = i)
      }, error = function(e) {
        cat("Error processing ensemble", i, model, scenario, year, ":", e$message, "\n")
        NULL
      })
    }
  }
}

# Aggregate decadal outputs
aggregate_decadal_outputs()

# Process ensemble outputs
process_ensemble_outputs()

# Generate documentation
generate_documentation()

cat("\n===== PROCESSING COMPLETE =====\n")
cat("Output directory:", output_dir, "\n")
cat("Total processing time:", format(Sys.time() - global_start_time), "\n")
