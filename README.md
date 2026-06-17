# GLEXIS

GLEXIS is a global, multi-model archive of monthly climate-extreme-day counts and conditional magnitudes derived from daily bias-adjusted ISIMIP3b forcing. The release covers historical simulations (1981-2014) and SSP1-2.6, SSP3-7.0, and SSP5-8.5 projections (2015-2100) for five CMIP6 models.

## Indicators

- Cold-extreme days: daily minimum temperature below the spatial p01 threshold.
- Hot days: daily maximum temperature above the spatial p95 threshold.
- Heavy-rain days: precipitation above the spatial p95 threshold and at least 1 mm d-1.
- High-wind days: near-surface mean wind above the spatial p95 threshold and at least 8 m s-1.
- Standardized water-balance lower-tail days: annual within-cell standardized precipitation-minus-PET values below a deposited p05 layer from the legacy baseline water-balance/SPEI-3 threshold workflow.

The water-balance indicator is not canonical SPEI or a general drought metric. Hot days are not run-length heatwaves, heavy-rain days are not floods, and high-wind days are not storm objects.

## Repository contents

- `R/production/`: historical and future GLEXIS production scripts.
- `R/thresholds/`: archived threshold-generation scripts.
- `release_metadata/`: implementation definitions, dictionaries, audits, manifest, and checksums.
- `data/README.md`: description of the separately deposited native raster and tabular archive.

## Reproducibility

The authoritative implementation definitions are in `release_metadata/implementation_definition.csv`. Public labels and legacy source identifiers are mapped in `release_metadata/hazard_variant_definition.csv`. Targeted single-cell corrections and their complete scope are documented in `release_metadata/targeted_corrections.csv`.


## Data availability

The raster and tabular data products are deposited separately because they are too large for a source-code repository. See `data/README.md` and `release_metadata/release_file_manifest.csv` for the data-package structure.


## License

Code is released under the MIT License. Data licensing is stated separately in the formal data repository record.
