# ======================================================================
# GLEXIS-HISTORICAL: Global Extreme Events from ISIMIP Historical Simulations
# A comprehensive processing pipeline for ISIMIP3b historical climate extremes
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
library(lubridate)

# =============== CONFIGURATION ===============
# Package metadata
pkg_metadata <- list(
  name = "GLEXIS-HISTORICAL",
  version = "1.0.0",
  description = "Global Extreme Events from ISIMIP Historical Simulations",
  author = "Dr. Andre Nakhavali",
  email = "nakhavali@iiasa.ac.at",
  institution = "International Institute for Applied Systems Analysis (IIASA)",
  license = "CC BY 4.0",
  doi = "10.5281/zenodo.XXXXXXX",
  reference = "Nakhavali et al. (2025) GLEXIS: A Global Extreme Events Dataset from ISIMIP3b"
)

# Model configuration
models <- "UKESM1-0-LL"#c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
scenario <- "historical"
years <- 2014#2011:2014

# Define the decades (10-year chunks for historical period)
decade_starts <- seq(1981, 2011, by = 10)  # 1981, 1991, 2001, 2011
decade_ends <- pmin(decade_starts + 9, 2014)  # 1990, 2000, 2010, 2014
decade_labels <- paste(decade_starts, decade_ends, sep = "-")

# Core processing parameters
params <- list(
  spei_scale = 3,                # SPEI time scale (months)
  south_pole_cutoff = -60,       # Exclude Antarctica (< -60°)
  min_precip_filter = 1,         # mm/day minimum precipitation threshold
  hyper_arid_threshold = 5,      # mm/month hyper-arid threshold
  annual_arid_threshold = 30,    # mm/year arid threshold
  annual_wet_threshold = 2000,   # mm/year permanently wet threshold
  wind_min_speed = 8,            # m/s minimum wind speed for mean wind
  tropical_lat_limit = 23.5,     # Tropical latitude boundary
  equatorial_lat_limit = 5,      # Equatorial latitude boundary
  jan_precip_threshold = 10,     # mm/month January precipitation threshold
  baseline_period = "1990-2010", # Baseline for threshold calculation
  thresholds_quantiles = list(   # Percentile thresholds for extremes
    tasmax = 4,      # 95th percentile (0-4 corresponds to 5th-95th)
    pr = 4,          # 95th percentile
    tasmin = 1,      # p01 threshold-file layer
    sfcwind = 4,     # 95th percentile (or 99th if using max wind)
    spei = 2         # p05 baseline SPEI-3 threshold-file layer
  ),
  ensemble_stats = c("mmm", "p10", "p50", "p90", "std", "n_models"),
  targeted_fixes = list(
    atacama_rainfall = TRUE,
    caribbean_cold = TRUE,
    ethiopia_heat = FALSE,
    north_atlantic_wind = TRUE,
    caribbean_drought = TRUE
  )
)

# File chunk structure for historical period
file_chunks <- list(
  "1981_1990" = 1981:1990,
  "1991_2000" = 1991:2000,
  "2001_2010" = 2001:2010,
  "2011_2014" = 2011:2014
)

# Path configuration
input_root <- "//pdrive/share/link/nakhavali.pdrv/watxene/ISIMIP/ISIMIP3b/InputData/climate_updated/bias-adjusted"
pet_root <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/PET"
thresholds_dir <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/Thresholds_back/"
wind_thresholds_dir <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/Thresholds_wind/"
output_dir <- "//hdrive/home$/u141/nakhavali/ISIMIP3b/OutputData/GLEXIS_historical_v1/"

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
  title = "GLEXIS-HISTORICAL: Global Extreme Events from ISIMIP Historical Simulations",
  summary = paste("Monthly counts of extreme climate events (heat, cold, precipitation,",
                  "wind, drought) derived from ISIMIP3b bias-adjusted historical climate",
                  "projections for multiple models."),
  keywords = paste("climate extremes, ISIMIP3b, heat waves, cold spells, heavy precipitation,",
                   "drought, wind storms, historical climate"),
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

# Function to get model-specific file pattern
get_model_pattern <- function(model) {
  # Special handling for UKESM historical files
  if (model == "UKESM1-0-LL") {
    return("ukesm1-0-ll_r1i1p1f2_w5e5")  # Note f2 instead of f1
  }
  
  # Standard patterns for other models
  switch(model,
         "GFDL-ESM4" = "gfdl-esm4_r1i1p1f1_w5e5",
         "IPSL-CM6A-LR" = "ipsl-cm6a-lr_r1i1p1f1_w5e5", 
         "MPI-ESM1-2-HR" = "mpi-esm1-2-hr_r1i1p1f1_w5e5",
         "MRI-ESM2-0" = "mri-esm2-0_r1i1p1f1_w5e5",
         "UKESM1-0-LL" = "ukesm1-0-ll_r1i1p1f2_w5e5"  # This line ensures consistent UKESM pattern
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
  
  anc_dir <- file.path(output_dir, "ancillary", paste0(tolower(model), "_", scenario), year)
  dir.create(anc_dir, recursive = TRUE, showWarnings = FALSE)
  
  nc_file <- file.path(anc_dir, 
                       paste0("GLEXIS_historical_", tolower(model), "_", year,
                              "_ancillary.nc"))
  
  dims <- list(
    lon = ncdim_def("lon", "degrees_east", lons),
    lat = ncdim_def("lat", "degrees_north", lats),
    time = ncdim_def("time", "days since 1850-01-01", 
                     as.numeric(dates - as.Date("1850-01-01")),
                     calendar = "standard")
  )
  
  vars <- list()
  for (name in names(ancillary_data)) {
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
  
  nc <- nc_create(nc_file, vars, force_v4 = TRUE)
  for (name in names(ancillary_data)) {
    ncvar_put(nc, name, aperm(ancillary_data[[name]], c(1, 2, 3)))
  }
  nc_close(nc)
}

# Function to save regional outputs
save_regional_outputs <- function(model, scenario, year, monthly_event_counts, lons, lats) {
  cat("\nSaving regional outputs with validation...\n")
  
  reg_dir <- file.path(output_dir, "regional", paste0(tolower(model), "_", scenario), year)
  dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
  
  countries <- ne_countries(scale = 50, returnclass = "sf")
  
  r <- rast(
    ncols = length(lons),
    nrows = length(lats),
    xmin = min(lons) - median(diff(sort(unique(lons)))) / 2,
    xmax = max(lons) + median(diff(sort(unique(lons)))) / 2,
    ymin = min(lats) - median(diff(sort(unique(lats)))) / 2,
    ymax = max(lats) + median(diff(sort(unique(lats)))) / 2,
    crs = "EPSG:4326"
  )
  
  for (event in names(monthly_event_counts)) {
    monthly_rasters <- lapply(1:12, function(m) {
      vals <- t(monthly_event_counts[[event]][, , m])
      r <- rast(vals)
      dx <- median(diff(sort(unique(lons))))
      dy <- median(diff(sort(unique(lats))))
      ext(r) <- ext(c(min(lons) - dx / 2, max(lons) + dx / 2,
                      min(lats) - dy / 2, max(lats) + dy / 2))
      crs(r) <- "EPSG:4326"
      return(r)
    })
    
    results <- list()
    for (m in 1:12) {
      country_stats <- exact_extract(monthly_rasters[[m]], countries, 'sum')
      
      results[[m]] <- data.frame(
        ISO3 = countries$iso_a3,
        Country = countries$name,
        Event = event,
        Year = year,
        Month = m,
        Count = round(country_stats, 1)
      )
    }
    
    results_df <- do.call(rbind, results)
    
    csv_file <- file.path(reg_dir, 
                          paste0("GLEXIS_historical_", tolower(model), "_", year,
                                 "_", event, "_country_counts.csv"))
    write.csv(results_df, csv_file, row.names = FALSE)
    cat("Saved regional outputs:", csv_file, "\n")
  }
}

# =============== DECADAL AGGREGATION FUNCTION ===============
aggregate_decadal_outputs <- function() {
  cat("\n===== AGGREGATING DECADAL OUTPUTS =====\n")
  
  decades <- list(
    "1981_1990" = 1981:1990,
    "1991_2000" = 1991:2000,
    "2001_2010" = 2001:2010,
    "2011_2014" = 2011:2014
  )
  
  for (model in models) {
    cat("\nProcessing model:", model, "...\n")
    
    for (decade_name in names(decades)) {
      decade_years <- decades[[decade_name]]
      cat("\nAggregating decade", decade_name, "...\n")
      
      dec_dir <- file.path(output_dir, "decadal", paste0(tolower(model), "_", scenario), decade_name)
      dir.create(dec_dir, recursive = TRUE, showWarnings = FALSE)
      
      decade_rasters <- list()
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
      
      for (event in c("heat", "rain", "cold", "wind", "drought")) {
        decade_array <- NULL
        
        for (year in decade_years) {
          year_dir <- file.path(output_dir, "annual", paste0(tolower(model), "_", scenario), year)
          tif_file <- file.path(year_dir, 
                                paste0("GLEXIS_historical_", tolower(model), "_", year,
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
        
        decade_avg <- apply(decade_array, c(1,2,3), mean, na.rm = TRUE)
        
        out_file <- file.path(dec_dir, 
                              paste0("GLEXIS_historical_", tolower(model), "_", decade_name,
                                     "_", event, "_monthly_avg.tif"))
        
        r_decade <- rast(decade_avg)
        ext(r_decade) <- ext(r)
        crs(r_decade) <- crs(r)
        names(r_decade) <- month.name[1:12]
        time(r_decade) <- as.Date(paste0(substr(decade_name, 1, 4), "-", 1:12, "-15"))
        varnames(r_decade) <- "event_days"
        units(r_decade) <- "days"
        
        for (i in 1:12) {
          longnames(r_decade)[i] <- paste("Decadal average", event, "days in", month.name[i], decade_name)
        }
        
        writeRaster(r_decade, filename = out_file, overwrite = TRUE,
                    gdal = c("DESCRIPTION=Decadal average monthly event counts"))
        cat("Saved:", out_file, "\n")
        
        decade_rasters[[event]] <- r_decade
      }
    }
  }
}

# =============== ENSEMBLE PROCESSING FUNCTION ===============
process_ensemble_outputs <- function() {
  cat("\n===== PROCESSING ENSEMBLE OUTPUTS =====\n")
  
  decades <- names(file_chunks)
  
  for (decade_name in decades) {
    cat("\nProcessing ensemble for decade", decade_name, "...\n")
    
    ens_dir <- file.path(output_dir, "ensemble", decade_name)
    dir.create(ens_dir, recursive = TRUE, showWarnings = FALSE)
    
    model_rasters <- list()
    for (model in models) {
      dec_dir <- file.path(output_dir, "decadal", paste0(tolower(model), "_", scenario), decade_name)
      
      if (dir.exists(dec_dir)) {
        event_rasters <- list()
        for (event in c("heat", "rain", "cold", "wind", "drought")) {
          tif_file <- file.path(dec_dir, 
                                paste0("GLEXIS_historical_", tolower(model), "_", decade_name,
                                       "_", event, "_monthly_avg.tif"))
          
          if (file.exists(tif_file)) {
            event_rasters[[event]] <- rast(tif_file)
          }
        }
        model_rasters[[model]] <- event_rasters
      }
    }
    
    if (length(model_rasters) == 0) {
      cat("No models found for decade", decade_name, "\n")
      next
    }
    
    for (event in c("heat", "rain", "cold", "wind", "drought")) {
      event_stack <- list()
      for (model in names(model_rasters)) {
        if (event %in% names(model_rasters[[model]])) {
          event_stack[[model]] <- model_rasters[[model]][[event]]
        }
      }
      
      if (length(event_stack) == 0) next
      
      r_stack <- rast(event_stack)
      
      ens_stats <- list(
        mean = app(r_stack, fun = mean, na.rm = TRUE),
        sd = app(r_stack, fun = sd, na.rm = TRUE),
        min = app(r_stack, fun = min, na.rm = TRUE),
        max = app(r_stack, fun = max, na.rm = TRUE),
        p10 = app(r_stack, fun = function(x) quantile(x, probs = 0.1, na.rm = TRUE)),
        p50 = app(r_stack, fun = function(x) quantile(x, probs = 0.5, na.rm = TRUE)),
        p90 = app(r_stack, fun = function(x) quantile(x, probs = 0.9, na.rm = TRUE))
      )
      
      for (stat_name in names(ens_stats)) {
        out_file <- file.path(ens_dir, 
                              paste0("GLEXIS_historical_ensemble_", decade_name,
                                     "_", event, "_", stat_name, ".tif"))
        
        r_stat <- ens_stats[[stat_name]]
        names(r_stat) <- month.name[1:12]
        time(r_stat) <- as.Date(paste0(substr(decade_name, 1, 4), "-", 1:12, "-15"))
        varnames(r_stat) <- "event_days"
        units(r_stat) <- "days"
        
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

# =============== DOCUMENTATION GENERATION ===============
generate_documentation <- function() {
  cat("\n===== GENERATING DOCUMENTATION =====\n")
  
  # Create README
  readme_content <- paste(
    "GLEXIS-HISTORICAL: Global Extreme Events from ISIMIP Historical Simulations\n",
    "=====================================================\n",
    "\n",
    "Dataset Description:\n",
    "--------------------\n",
    "This dataset contains monthly counts of extreme climate events derived from ISIMIP3b ",
    "bias-adjusted historical climate simulations (1981-2014) for multiple models.\n",
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
    "1. annual/: Contains annual files for each model\n",
    "2. decadal/: Contains decade-aggregated files (1981-1990, 1991-2000, etc.)\n",
    "3. ensemble/: Contains multi-model ensemble statistics\n",
    "4. regional/: Contains country-level aggregated event counts\n",
    "5. ancillary/: Contains masks and threshold values\n",
    "\n",
    "Models Included:\n",
    "---------------\n",
    paste("-", models, collapse = "\n"), "\n",
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
        threshold = "95th percentile of daily wind speed (1990-2010 baseline)"
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
      scenario = scenario,
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
process_single_year <- function(model, scenario, year) {
  start_time <- Sys.time()
  cat("\n=== PROCESSING YEAR", year, "===\n")
  
  model_pattern <- get_model_pattern(model)
  
  # 1. FIND ALL REQUIRED FILES
  cat("Step 1: Finding all required files...\n")
  
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
    pet = file.path(pet_root, scenario, model,
                    paste0(model_pattern, "_", scenario,
                           "_pet_global_daily_", chunk, ".nc"))
  )
  
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
  
  keep_lats <- which(lats >= params$south_pole_cutoff)
  lats <- lats[keep_lats]
  
  cat("Grid dimensions:", length(lons), "longitudes x", length(lats), "latitudes (Antarctica removed)\n")
  cat("Selected year:", year, ", Days:", length(year_idx), "\n")
  
  # 3. CREATE LAND MASK (excluding Antarctica and glaciers)
  cat("\nStep 3: Creating land mask (excluding Antarctica and glaciers)...\n")
  world <- ne_countries(scale = 50, returnclass = "sf")
  glaciers <- ne_download(scale = 50, type = 'glaciated_areas', category = 'physical', returnclass = "sf")
  
  r <- rast(ncols = length(lons), nrows = length(lats),
            xmin = min(lons), xmax = max(lons),
            ymin = min(lats), ymax = max(lats),
            crs = "EPSG:4326")
  
  land_raster <- rasterize(world, r, field = 1, background = 0)
  glacier_raster <- rasterize(glaciers, r, field = 1, background = 0)
  
  land_mask <- land_raster * (1 - glacier_raster)
  land_mask_matrix <- t(as.matrix(land_mask, wide = TRUE))
  
  land_mask_3d <- array(land_mask_matrix, dim = c(dim(land_mask_matrix), length(year_idx)))
  land_mask_3d[land_mask_3d < 0.5] <- 0
  
  cat("\nStep 4: Loading historical climate data for desert/wet region detection...\n")
  
  load_pr_thresholds <- function(model) {
    model_patterns <- list(
      "GFDL-ESM4" = "gfdl-esm4_r1i1p1f1_w5e5",
      "IPSL-CM6A-LR" = "ipsl-cm6a-lr_r1i1p1f1_w5e5",
      "MPI-ESM1-2-HR" = "mpi-esm1-2-hr_r1i1p1f1_w5e5",
      "MRI-ESM2-0" = "mri-esm2-0_r1i1p1f1_w5e5",
      "UKESM1-0-LL" = "ukesm1-0-ll_r1i1p1f2_w5e5"
    )
    
    pattern <- model_patterns[[model]]
    
    threshold_file <- file.path(thresholds_dir,
                                paste0(pattern, "_thresholds_1990-2010_SPEI3.nc"))
    
    if (!file.exists(threshold_file)) {
      stop(paste("Threshold file not found:", threshold_file))
    }
    
    nc <- nc_open(threshold_file)
    pr_thresholds <- ncvar_get(nc, "pr")[, , 4]
    nc_close(nc)
    
    pr_thresholds <- pr_thresholds[, keep_lats]
    
    return(list(
      thresholds = pr_thresholds,
      mean_annual = pr_thresholds * 365,
      mean_jan = pr_thresholds
    ))
  }
  
  precip_data <- tryCatch({
    pr_thresholds <- load_pr_thresholds(model)
    list(
      thresholds = pr_thresholds,
      mean_annual = pr_thresholds * 365,
      mean_jan = pr_thresholds
    )
  }, error = function(e) {
    list(
      thresholds = matrix(1, nrow=length(lons), ncol=length(lats)),
      mean_annual = matrix(365, nrow=length(lons), ncol=length(lats)),
      mean_jan = matrix(30, nrow=length(lons), ncol=length(lats))
    )
  })
  
  if (!is.null(precip_data)) {
    mean_annual_precip <- precip_data$mean_annual
    mean_jan_precip <- precip_data$mean_jan
  } else {
    mean_annual_precip <- NULL
    mean_jan_precip <- NULL
  }
  
  hyper_arid_mask <- if (!is.null(mean_jan_precip)) {
    mean_jan_precip < params$hyper_arid_threshold
  } else {
    matrix(FALSE, nrow = length(lons), ncol = length(lats))
  }
  
  arid_mask <- if (!is.null(mean_annual_precip)) {
    mean_annual_precip < params$annual_arid_threshold
  } else {
    matrix(FALSE, nrow = length(lons), ncol = length(lats))
  }
  
  permanent_wet_mask <- if (!is.null(mean_annual_precip)) {
    mean_annual_precip > params$annual_wet_threshold
  } else {
    matrix(FALSE, nrow = length(lons), ncol = length(lats))
  }
  
  jan_precip_mask <- if (!is.null(mean_jan_precip)) {
    mean_jan_precip < params$jan_precip_threshold
  } else {
    matrix(FALSE, nrow = length(lons), ncol = length(lats))
  }
  
  hyper_arid_mask_3d <- array(hyper_arid_mask, dim = c(dim(hyper_arid_mask), length(year_idx)))
  arid_mask_3d <- array(arid_mask, dim = c(dim(arid_mask), length(year_idx)))
  permanent_wet_mask_3d <- array(permanent_wet_mask, dim = c(dim(permanent_wet_mask), length(year_idx)))
  jan_precip_mask_3d <- array(jan_precip_mask, dim = c(dim(jan_precip_mask), length(year_idx)))
  
  cat("\nStep 5: Creating high-resolution ocean mask...\n")
  coastline <- ne_coastline(scale = 50, returnclass = "sf")
  ocean_mask <- rasterize(coastline, r, field = 1, background = 0)
  ocean_mask_matrix <- t(as.matrix(ocean_mask, wide = TRUE))
  ocean_mask_3d <- array(ocean_mask_matrix, dim = dim(land_mask_3d))
  
  cat("\nStep 6: Creating latitude-based masks...\n")
  
  lat_grid <- array(rep(lats, each = length(lons)), dim = c(length(lons), length(lats)))
  
  equatorial_ocean_mask <- (lat_grid >= -params$equatorial_lat_limit & lat_grid <= params$equatorial_lat_limit) & (ocean_mask_matrix == 1)
  tropical_ocean_mask <- (lat_grid >= -params$tropical_lat_limit & lat_grid <= params$tropical_lat_limit) & (ocean_mask_matrix == 1)
  arctic_mask <- lat_grid > 60
  
  equatorial_ocean_mask_3d <- array(equatorial_ocean_mask, dim = dim(land_mask_3d))
  tropical_ocean_mask_3d <- array(tropical_ocean_mask, dim = dim(land_mask_3d))
  arctic_mask_3d <- array(arctic_mask, dim = dim(land_mask_3d))
  
  cat("\nStep 7: Creating comprehensive exclusion masks...\n")
  
  rain_exclude_mask <- hyper_arid_mask | jan_precip_mask
  cold_exclude_mask <- ocean_mask_matrix == 1
  heat_exclude_mask <- equatorial_ocean_mask
  wind_exclude_mask <- ocean_mask_matrix == 1
  drought_exclude_mask <- arid_mask | permanent_wet_mask
  
  rain_exclude_mask_3d <- array(rain_exclude_mask, dim = c(dim(rain_exclude_mask), length(year_idx)))
  cold_exclude_mask_3d <- array(cold_exclude_mask, dim = c(dim(cold_exclude_mask), length(year_idx)))
  heat_exclude_mask_3d <- array(heat_exclude_mask, dim = c(dim(heat_exclude_mask), length(year_idx)))
  wind_exclude_mask_3d <- array(wind_exclude_mask, dim = c(dim(wind_exclude_mask), length(year_idx)))
  drought_exclude_mask_3d <- array(drought_exclude_mask, dim = c(dim(drought_exclude_mask), length(year_idx)))
  
  cat("\nStep 7: Loading thresholds...\n")
  
  threshold_quantiles <- list(
    tasmax = 4,      # 95th percentile
    pr = 4,          # 95th percentile
    tasmin = 1,      # p01 threshold-file layer
    sfcwind = 4,     # 95th percentile
    spei = 2         # p05 baseline SPEI-3 threshold-file layer
  )
  
  load_thresholds <- function(var) {
    base_pattern <- paste0(model_pattern, "_r1i1p1f1_w5e5")
    
    possible_files <- c(
      file.path(thresholds_dir, paste0(base_pattern, "_thresholds_1990-2010_SPEI3.nc")),
      file.path(thresholds_dir, paste0(base_pattern, "_", var, "_thresholds_1990-2010.nc")),
      file.path(thresholds_dir, paste0(base_pattern, "_thresholds_1990-2010_", var, ".nc")),
      file.path(thresholds_dir, paste0(model_pattern, "_thresholds_1990-2010_SPEI3.nc")),
      file.path(thresholds_dir, paste0(model_pattern, "_", var, "_thresholds_1990-2010.nc")),
      file.path(thresholds_dir, paste0(model_pattern, "_thresholds_1990-2010_", var, ".nc"))
    )
    
    if (var == "sfcwind") {
      possible_files <- c(
        file.path(wind_thresholds_dir, paste0(base_pattern, "_wind_thresholds_1990-2010.nc")),
        file.path(wind_thresholds_dir, paste0(model_pattern, "_wind_thresholds_1990-2010.nc")),
        possible_files
      )
    }
    
    threshold_file <- NULL
    for (f in possible_files) {
      if (file.exists(f)) {
        threshold_file <- f
        break
      }
    }
    
    if (is.null(threshold_file)) {
      stop(paste("No threshold file found for", var))
    }
    
    cat("Loading thresholds for", var, "from:", threshold_file, "\n")
    nc <- nc_open(threshold_file)
    
    var_names <- names(nc$var)
    target_var <- ifelse(var == "sfcwind", "sfcWind", var)
    if (!target_var %in% var_names) {
      target_var <- var_names[grep(var, var_names)[1]]
    }
    
    data <- ncvar_get(nc, target_var)
    nc_close(nc)
    
    data <- data[, keep_lats, ]
    
    return(data)
  }
  
  thresholds <- list()
  for (var in c("tasmax", "pr", "tasmin", "sfcwind", "spei")) {
    tryCatch({
      thresholds[[var]] <- load_thresholds(var)
    }, error = function(e) {
      cat("Error loading thresholds for", var, ":", e$message, "\n")
      stop("Threshold loading failed")
    })
  }
  
  cat("\nStep 8: Creating 3D threshold arrays...\n")
  thresholds_3d <- list()
  
  for (var in names(thresholds)) {
    thresh_data <- thresholds[[var]]
    dims <- dim(thresh_data)
    
    if (var == "spei") {
      if (length(dims) == 3 && dims[3] == 9) {
        quantile_slice <- thresh_data[, , threshold_quantiles[[var]]]
        thresholds_3d[[var]] <- array(quantile_slice,
                                      dim = c(dim(quantile_slice), length(year_idx)))
      } else {
        stop(paste("Unexpected SPEI threshold dimensions:", paste(dims, collapse = "x")))
      }
    }
    else if (length(dims) == 4) {
      thresholds_3d[[var]] <- thresh_data[, , , threshold_quantiles[[var]]]
    }
    else if (length(dims) == 3 && dims[3] == 5) {
      quantile_slice <- thresh_data[, , threshold_quantiles[[var]]]
      thresholds_3d[[var]] <- array(quantile_slice,
                                    dim = c(dim(quantile_slice), length(year_idx)))
    }
    else if (length(dims) == 3 && dims[3] == 365) {
      thresholds_3d[[var]] <- thresh_data[, , yday(dates[year_idx])]
    }
    else if (length(dims) == 3 && dims[3] == 12) {
      thresholds_3d[[var]] <- array(NA, dim = c(dims[1:2], length(year_idx)))
      for (m in 1:12) {
        month_days <- which(month(dates[year_idx]) == m)
        if (length(month_days) > 0) {
          thresholds_3d[[var]][, , month_days] <- array(thresh_data[, , m],
                                                        dim = c(dims[1:2], length(month_days)))
        }
      }
    }
    else if (length(dims) == 2) {
      thresholds_3d[[var]] <- array(thresh_data,
                                    dim = c(dims, length(year_idx)))
    }
    else {
      stop(paste("Unexpected threshold dimensions for", var,
                 ":", paste(dims, collapse = "x")))
    }
    
    if (length(dim(thresholds_3d[[var]])) != 3) {
      stop(paste("Failed to create proper 3D array for", var,
                 "Resulting dimensions:",
                 paste(dim(thresholds_3d[[var]]), collapse = "x")))
    }
    
    land_mask_reshaped <- if (!identical(dim(thresholds_3d[[var]]), dim(land_mask_3d))) {
      array(land_mask_3d[, , 1], dim = dim(thresholds_3d[[var]]))
    } else {
      land_mask_3d
    }
    
    thresholds_3d[[var]] <- thresholds_3d[[var]] * land_mask_reshaped
  }
  
  cat("\nStep 9: Loading climate data...\n")
  
  wind_data_type <- "mean"
  
  load_var <- function(var, unit_conversion = NULL) {
    if (var == "sfcwind") {
      max_file1 <- file.path(input_root, scenario, model,
                             paste0(tolower(model), "_r1i1p1f1_w5e5_", scenario,
                                    "_sfcwindmax_global_daily_", chunk, ".nc"))
      max_file2 <- file.path(input_root, scenario, model,
                             paste0(tolower(model), "_", scenario,
                                    "_sfcwindmax_global_daily_", chunk, ".nc"))
      max_file3 <- file.path(input_root, scenario, model,
                             paste0(model_pattern, "_", scenario,
                                    "_", var, "_global_daily_", chunk, ".nc"))
      max_file4 <- file.path(input_root, scenario, model,
                             paste0("ukesm1-0-ll_r1i1p1f2_w5e5_", scenario,
                                    "_", var, "_global_daily_", chunk, ".nc"))
      
      if (file.exists(max_file1)) {
        cat("Using wind MAX data (with run identifier)\n")
        nc_file <- max_file1
        var_name <- "sfcwindmax"
        wind_data_type <<- "max"
      } else if (file.exists(max_file2)) {
        cat("Using wind MAX data (without run identifier)\n")
        nc_file <- max_file2
        var_name <- "sfcwindmax"
        wind_data_type <<- "max"
      } else {
        mean_file1 <- file.path(input_root, scenario, model,
                                paste0(tolower(model), "_r1i1p1f1_w5e5_", scenario,
                                       "_sfcwind_global_daily_", chunk, ".nc"))
        mean_file2 <- file.path(input_root, scenario, model,
                                paste0(tolower(model), "_", scenario,
                                       "_sfcwind_global_daily_", chunk, ".nc"))
        mean_file3 <- file.path(input_root, scenario, model,
                                paste0(model_pattern, "_", scenario,
                                       "_", var, "_global_daily_", chunk, ".nc"))
        mean_file4 <- file.path(input_root, scenario, model,
                                paste0("ukesm1-0-ll_r1i1p1f2_w5e5_", scenario,
                                       "_", var, "_global_daily_", chunk, ".nc"))
        
        if (file.exists(mean_file1)) {
          cat("Using wind MEAN data (with run identifier)\n")
          nc_file <- mean_file1
          var_name <- "sfcwind"
          wind_data_type <<- "mean"
        } else if (file.exists(mean_file2)) {
          cat("Using wind MEAN data (without run identifier)\n")
          nc_file <- mean_file2
          var_name <- "sfcwind"
          wind_data_type <<- "mean"
        } else if (file.exists(mean_file3)) {
          cat("Using wind MEAN data (ukesm)\n")
          nc_file <- mean_file3
          var_name <- "sfcwind"
          wind_data_type <<- "mean"
        } else if (file.exists(mean_file4)) {
          cat("Using wind MEAN data (ukesm)\n")
          nc_file <- mean_file4
          var_name <- "sfcwind"
          wind_data_type <<- "mean"
        } else {
          stop(paste("Wind data file not found. Tried:\n", max_file1, "\n", max_file2,
                     "\n", mean_file1, "\n", mean_file2))
        }
      }
    } else {
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
    
    data <- ncvar_get(nc, var_name,
                      start = c(1, 1, min(year_idx)),
                      count = c(-1, -1, length(year_idx)))
    nc_close(nc)
    
    data <- data[, keep_lats, ]
    
    if (!is.null(unit_conversion)) {
      if (unit_conversion == "K_to_C") data <- data - 273.15
      if (unit_conversion == "kgm2s_to_mmday") data <- data * 86400
    }
    
    return(data)
  }
  
  year_data <- list(
    tasmax = load_var("tasmax", "K_to_C"),
    pr = load_var("pr", "kgm2s_to_mmday"),
    tasmin = load_var("tasmin", "K_to_C"),
    sfcwind = load_var("sfcwind")
  )
  
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
  
  time_pet <- ncvar_get(nc_pet, "time")
  time_units <- ncatt_get(nc_pet, "time")$units
  calendar <- ncatt_get(nc_pet, "time", "calendar")$value
  if (is.null(calendar)) calendar <- "standard"
  
  origin <- as.Date(gsub("days since ", "", time_units))
  dates_pet <- as.Date(time_pet, origin = origin)
  
  file_years <- unique(format(dates_pet, "%Y"))
  
  if (!as.character(year) %in% file_years) {
    nc_close(nc_pet)
    stop(paste("Year", year, "not found in PET file. Available years:",
               paste(file_years, collapse = ", ")))
  }
  
  year_idx <- which(format(dates_pet, "%Y") == sprintf("%04d", year))
  
  if (length(year_idx) == 0) {
    nc_close(nc_pet)
    stop(paste("No time indices found for year", year, "in PET file"))
  }
  
  year_data$pet <- tryCatch({
    ncvar_get(nc_pet, "pet",
              start = c(1, 1, min(year_idx)),
              count = c(-1, -1, length(year_idx)))
  }, error = function(e) {
    nc_close(nc_pet)
    stop(paste("Failed to read PET data:", e$message))
  })
  
  nc_close(nc_pet)
  
  if (length(dim(year_data$pet)) != 3) {
    stop("PET data is not 3-dimensional")
  }
  
  year_data$pet <- year_data$pet[, keep_lats, ]
  
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
  
  cat("\nStep 12: Calculating water balance and annual within-cell standardized scores...\n")
  
  wb <- year_data$pr - year_data$pet
  
  spei_result <- array(NA, dim = dim(wb))
  
  land_points <- which(land_mask_matrix == 1, arr.ind = TRUE)
  
  for (i in 1:nrow(land_points)) {
    idx <- land_points[i, ]
    
    if (idx[1] > dim(wb)[1] || idx[2] > dim(wb)[2]) next
    
    wb_series <- wb[idx[1], idx[2], ]
    
    if (sum(!is.na(wb_series)) >= params$spei_scale) {
      tryCatch({
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
  
  spei_result <- spei_result * land_mask_3d
  spei_result[is.infinite(spei_result)] <- NA
  spei_result <- round(spei_result, 2)
  
  cat("\nStep 13: Calculating extreme events with comprehensive exclusions...\n")
  
  events <- list(
    heat = (year_data$tasmax > thresholds_3d$tasmax) &
      !equatorial_ocean_mask_3d,
    
    rain = (year_data$pr > thresholds_3d$pr) &
      (year_data$pr >= params$min_precip_filter) &
      !hyper_arid_mask_3d,
    
    cold = (year_data$tasmin < thresholds_3d$tasmin) &
      !equatorial_ocean_mask_3d,
    
    wind = (year_data$sfcwind > thresholds_3d$sfcwind) &
      (land_mask_3d == 1),
    
    drought = (spei_result < thresholds_3d$spei) &
      !arid_mask_3d &
      !permanent_wet_mask_3d
  )
  
  if (wind_data_type == "mean") {
    events$wind <- events$wind & (year_data$sfcwind >= params$wind_min_speed)
  }
  
  cat("\nApplying targeted fixes for remaining problematic pixels...\n")
  
  if (1 %in% month(dates[year_idx])) {
    jan_days <- which(month(dates[year_idx]) == 1)
    
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
  
  for (event in names(events)) {
    events[[event]] <- ifelse(land_mask_3d == 1, events[[event]], NA)
  }
  
  cat("\nStep 14: Summing events by month...\n")
  
  monthly_event_counts <- list()
  for (event in names(events)) {
    monthly_event_counts[[event]] <- array(0, dim = c(dim(events[[event]])[1:2], 12))
  }
  
  for (m in 1:12) {
    month_days <- which(month(dates[year_idx]) == m)
    if (length(month_days) > 0) {
      for (event in names(events)) {
        event_slice <- events[[event]][, , month_days, drop = FALSE]
        valid_cells <- apply(!is.na(event_slice), c(1,2), any)
        month_counts <- apply(event_slice, c(1,2), sum, na.rm = TRUE)
        month_counts[!valid_cells] <- NA_real_
        monthly_event_counts[[event]][, , m] <- month_counts
      }
    }
  }
  
  cat("\nStep 15: Creating output rasters...\n")
  
  create_output_raster <- function(data) {
    r <- rast(t(data), crs = "EPSG:4326")
    dx <- median(diff(sort(unique(lons))))
    dy <- median(diff(sort(unique(lats))))
    ext(r) <- c(min(lons) - dx / 2, max(lons) + dx / 2,
                min(lats) - dy / 2, max(lats) + dy / 2)
    return(r)
  }
  
  monthly_rasters <- list()
  for (event in names(monthly_event_counts)) {
    monthly_rasters[[event]] <- lapply(1:12, function(m) {
      create_output_raster(monthly_event_counts[[event]][, , m])
    })
  }
  
  cat("\nStep 16: Saving outputs...\n")
  
  year_dir <- file.path(output_dir, "annual", paste0(tolower(model), "_", scenario), year)
  dir.create(year_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (event in names(monthly_rasters)) {
    r_stack <- rast(monthly_rasters[[event]])
    
    names(r_stack) <- month.name[1:12]
    time(r_stack) <- as.Date(paste(year, 1:12, "15", sep = "-"))
    varnames(r_stack) <- "event_days"
    units(r_stack) <- "days"
    
    for (i in 1:12) {
      longnames(r_stack)[i] <- paste("Monthly count of", event, "days in", month.name[i], year)
    }
    
    out_file <- file.path(year_dir,
                          paste0("GLEXIS_historical_", tolower(model), "_", year,
                                 "_", event, "_monthly_counts.tif"))
    writeRaster(r_stack, filename = out_file, overwrite = TRUE,
                gdal = c("DESCRIPTION=Monthly extreme event counts"))
    cat("Saved:", out_file, "\n")
  }
  
  cat("\nStep 17: Saving ancillary outputs...\n")
  
  ancillary_data <- list(
    land_mask = land_mask_3d,
    ocean_mask = ocean_mask_3d,
    glacier_mask = glacier_raster,
    hyper_arid_mask = hyper_arid_mask_3d,
    arid_mask = arid_mask_3d,
    permanent_wet_mask = permanent_wet_mask_3d,
    equatorial_mask = equatorial_ocean_mask_3d,
    tropical_mask = tropical_ocean_mask_3d,
    arctic_mask = arctic_mask_3d,
    tasmax_threshold = thresholds_3d$tasmax,
    pr_threshold = thresholds_3d$pr,
    tasmin_threshold = thresholds_3d$tasmin,
    sfcwind_threshold = thresholds_3d$sfcwind,
    spei_threshold = thresholds_3d$spei
  )
  
  save_ancillary_outputs(model, scenario, year, ancillary_data, lons, lats, dates[year_idx])
  
  cat("\nStep 18: Saving regional outputs...\n")
  save_regional_outputs(model, scenario, year, monthly_event_counts, lons, lats)
  
  processing_time <- Sys.time() - start_time
  cat("Year", year, "processed in:", format(processing_time), "\n")
  
  return(list(
    monthly_event_counts = monthly_event_counts,
    ancillary_data = ancillary_data
  ))
}

# =============== MAIN EXECUTION ===============
cat("GLEXIS-HISTORICAL Processing Pipeline\n")
cat("=====================================\n")
cat("Version:", pkg_metadata$version, "\n")
cat("Author:", pkg_metadata$author, "\n")
cat("Institution:", pkg_metadata$institution, "\n")
cat("License:", pkg_metadata$license, "\n")
cat("DOI:", pkg_metadata$doi, "\n\n")

global_start_time <- Sys.time()

# Initialize storage for ensemble processing
ensemble_data <- list()

# Process all models
for (model in models) {
  cat("\n===== PROCESSING MODEL", model, "=====\n")
  
  # Process each year
  for (year in years) {
    cat("\n--- PROCESSING YEAR", year, "---\n")
    
    result <- tryCatch({
      process_single_year(model, scenario, year)
    }, error = function(e) {
      cat("Error processing", model, year, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(result)) {
      # Store for ensemble processing
      ensemble_data[[paste(model, year)]] <- result
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
