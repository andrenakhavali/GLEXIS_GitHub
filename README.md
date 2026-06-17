# GLEXIS

**GLEXIS** is the **Global Extreme event Indicators from ISIMIP Simulations** dataset and code archive. It provides a global, multi-model set of monthly climate-extreme indicators derived from daily bias-adjusted ISIMIP3b atmospheric forcing for 1981-2100.

The dataset is designed for comparative exposure, impact-model, and adaptation-analysis workflows that need spatially explicit, multi-hazard, multi-model climate-extreme information on a common historical baseline.

## Scope

GLEXIS covers:

- Historical simulations for 1981-2014.
- Future projections for 2015-2100 under SSP1-2.6, SSP3-7.0, and SSP5-8.5.
- Five CMIP6 models: GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, and UKESM1-0-LL.
- Valid land grid cells on the native 0.5 degree ISIMIP3b grid north of 60 degrees S.
- Monthly raster products and country-level tabular summaries for selected products.

Thresholds are calibrated to a fixed 1990-2010 historical baseline and then applied consistently to historical and future periods. This fixed-baseline design makes changes in event frequency interpretable relative to the same reference climatology.

## Indicators

GLEXIS includes five headline meteorological hazard indicators:

- **Hot days**: daily maximum near-surface air temperature above the local historical p95 threshold. This is a daily frequency indicator, not a heatwave spell metric.
- **Cold-extreme days**: daily minimum near-surface air temperature below the local historical p01 threshold. This is not a fixed-threshold frost-day index.
- **Heavy-rain days**: daily precipitation above the local historical p95 threshold and at least 1 mm d-1. This is a meteorological precipitation indicator, not a flood indicator.
- **High-wind days**: daily mean near-surface wind speed above the local historical p95 threshold and at least 8 m s-1. This uses daily mean wind speed, not observed gusts.
- **Standardized water-balance drought indicator**: annual precipitation-minus-PET water balance standardized against the 1990-2010 baseline and classified below the local p05 threshold. This is an annual water-balance deficit indicator, not canonical SPEI-3.

The drought indicator is annual in resolution. In the monthly raster archive, a drought-classified year is represented uniformly across all 12 months for archive consistency.

## Products

The full data archive is deposited separately from this source-code repository because the raster products are too large for GitHub.

The manuscript data record describes:

- **4,080 historical GeoTIFF files**: 5 models x 34 years x 24 files.
- **30,960 future GeoTIFF files**: 5 models x 3 SSPs x 86 years x 24 files.
- Twelve monthly bands per GeoTIFF, January through December.
- WGS84 spatial reference (EPSG:4326), 0.5 degree pixel size, and 720 x 300 raster dimensions.
- Event-count layers and conditional minimum, mean, and maximum magnitude layers.
- Historical country-level annual and monthly ensemble CSV products for 178 countries.
- Global area-weighted annual summaries across historical and future periods.

Country-level CSV column definitions are documented in `release_metadata/csv_data_dictionary.csv`.

## Repository Contents

- `R/production/`: production scripts for historical and future GLEXIS processing.
- `R/thresholds/`: archived threshold-generation scripts.
- `release_metadata/`: implementation definitions, data dictionaries, file manifests, correction logs, validation summaries, and software provenance.
- `data/README.md`: notes on the separately deposited data archive.
- `CITATION.cff`: software citation metadata for GitHub and Zenodo-style releases.
- `LICENSE`: MIT license for the code in this repository.

## Reproducibility

The authoritative implementation details are recorded in:

- `release_metadata/implementation_definition.csv`
- `release_metadata/hazard_variant_definition.csv`
- `release_metadata/targeted_corrections.csv`
- `release_metadata/calendar_harmonization_audit.csv`
- `release_metadata/software_provenance.csv`

Key reproducibility conventions from the manuscript:

- ISIMIP3b daily forcing is bias-adjusted against W5E5.
- Thresholds use a fixed 1990-2010 baseline.
- Gregorian leap years are harmonized by removing the nominal 29 February slot for consistent day-of-year indexing.
- UKESM1-0-LL keeps its native 360-day calendar.
- GeoTIFF count-layer NoData is encoded as `-9999`; conditional magnitude NoData is encoded as `NaN`.
- Six targeted single-cell/month corrections are documented in the release metadata.

## Usage Notes

GLEXIS is suitable for:

- Monthly multi-hazard forcing for impact models.
- Multi-model and multi-scenario comparison of climate-hazard frequency trajectories.
- Cross-hazard co-occurrence and compound-extreme screening.
- Historical trend characterization relative to the fixed 1990-2010 baseline.

Important limitations:

- GLEXIS indicators are meteorological hazard screening metrics, not impact or risk estimates.
- Free-running GCM simulations should not be compared year-by-year with specific observed disaster events.
- The high-wind indicator uses daily mean wind speed; wind-damage studies should treat it carefully.
- Wind and drought validation confidence is lower than for temperature and precipitation indicators.
- The ensemble mean may not represent any individual model well in regions with high inter-model spread.

## Data Availability

The main GLEXIS data-link browser is available at:

https://andrenakhavali.github.io/GLEXIS_data_links_public/

The raster GeoTIFF files, tabular CSV products, threshold NetCDF files, and full release metadata are deposited separately from this code repository. Add the final data DOI here once the formal data record is published.

See `data/README.md` and `release_metadata/release_file_manifest.csv` for the expected data-package structure.

## Code Availability

The production and threshold-generation scripts are released in this repository under the MIT License. The main production scripts are:

- `R/production/GLEXIS_historical.R`
- `R/production/GLEXIS_future_with_wind_gust.R`
- `R/production/extreme_thresholds_wind_with_gust.R`

## Citation

Please cite the GLEXIS data descriptor manuscript and the archived software release. Repository-level citation metadata are provided in `CITATION.cff`.

## License

Code is released under the MIT License. Data licensing is stated separately in the formal data repository record.
