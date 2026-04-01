# Causal Impact of London's Ultra Low Emission Zone (ULEZ) on Air Quality

Evaluating the effect of London's staggered ULEZ expansion on NO₂ and NOₓ concentrations using **Random Forest weather normalization** and **partially pooled Synthetic Control Methods** (Ben-Michael et al., 2021).

## Research Question

Did London's Ultra Low Emission Zone — implemented in three staggered phases (Central 2019, Inner 2021, Outer 2023) — cause a statistically significant reduction in roadside NO₂ and NOₓ concentrations, after removing the confounding effects of weather variability?

## Key Findings

| Pollutant | Average Treatment Effect (ATT) | Percentage Change | p-value (Placebo) |
|-----------|-------------------------------|-------------------|-------------------|
| NO₂       | −0.054 (log scale)            | ≈ −5.4%           | < 0.05            |
| NOₓ       | −0.064 (log scale)            | ≈ −6.4%           | < 0.05            |

Results are robust to spatial placebo tests, in-time placebo tests, and leave-one-out cross-validation.

## Methodology

### Pipeline Overview

```
Hourly pollution data (KCL / AURN / WAQN / NI APIs)
        │
        ▼
Missing value imputation (Kalman smoothing, maxgap = 48h)
        │
        ▼
Merge with weather data (NOAA + ERA5 reanalysis)
        │
        ▼
Random Forest weather normalization (rmweather)
   ├── 500 trees per site
   ├── Predictors: wind speed/direction, temperature, humidity,
   │   pressure, BLH, solar radiation, cloud cover, precipitation,
   │   and temporal features (hour, weekday, month, day_julian)
   └── OOB R² typically 0.77–0.88
        │
        ▼
Weekly aggregation of weather-normalized concentrations
        │
        ▼
Staggered Adoption Synthetic Control (augsynth::multisynth)
   ├── Log-transformed outcomes
   ├── Optimal nu selection via balance frontier
   ├── 52-week pre/post treatment window
   └── Unit fixed effects (fixedeff = TRUE)
        │
        ▼
Robustness checks
   ├── Spatial placebo tests (500 permutations)
   ├── In-time placebo tests
   └── Leave-one-out cross-validation
```

### Weather Normalization

Air pollution concentrations are heavily influenced by meteorological conditions. To isolate the policy effect, I use Random Forest models (via the `rmweather` package) to predict hourly pollution concentrations from weather and temporal variables, then generate counterfactual predictions by randomly shuffling weather variables 500 times. The resulting weather-normalized series removes meteorological noise while preserving the underlying emission trend.

### Staggered Synthetic Control

London's ULEZ was rolled out in three phases, creating a staggered adoption setting. I use the `augsynth` package's `multisynth` function, which constructs a weighted combination of control units (38 monitoring sites across UK cities outside London) to approximate the counterfactual trajectory for each treated unit. The partially pooled estimator balances between unit-specific and globally pooled synthetic controls, with the pooling parameter (nu) selected via a balance frontier heuristic.

## Data Sources

| Data | Source | Access Method |
|------|--------|---------------|
| London pollution (NO₂, NOₓ) | King's College London (KCL) | `openair::importKCL()` |
| UK pollution (control sites) | AURN / WAQN / NI networks | `openair::importAURN()`, `importWAQN()`, `importNI()` |
| Surface weather (London only) | NOAA Integrated Surface Database | `worldmet::importNOAA()` |
| Weather variables (all sites) | ERA5 reanalysis (ECMWF) | CDS API (Python script) |

Weather variables from ERA5 include: wind speed, wind direction, air temperature, relative humidity, surface pressure, total cloud cover, boundary layer height, surface solar radiation, and total precipitation.

**Note:** Raw data files are not included in this repository due to size. All data can be reproduced using the provided scripts — see [Reproducing the Analysis](#reproducing-the-analysis).

## Repository Structure

```
london-ulez-air-quality/
│
├── README.md
├── .gitignore
│
├── config/
│   └── sites_config.csv              # Site metadata and treatment assignment
│
├── scripts/
│   ├── 01_download_era5_weather.py   # Download ERA5 data via CDS API
│   ├── 02_prepare_weather.R          # Merge NOAA + ERA5, interpolate missing
│   ├── 03_import_and_normalize.R     # Import pollution → fill gaps → RF normalize
│   ├── 04_build_panel.R              # Consolidate weekly panel data
│   ├── 05_synthetic_control.R        # Nu selection + multisynth estimation
│   ├── 06_robustness_checks.R        # Placebo tests + LOOCV
│   └── utils/
│       ├── fill_missing.R            # Kalman smoothing imputation
│       ├── analyze_site.R            # Weather normalization core function
│       └── plotting_helpers.R        # Visualization functions
│
├── data/
│   ├── README.md                     # Data dictionary and download instructions
│   └── site_metadata/
│       └── london_valid_sites.csv
│
├── outputs/
│   ├── figures/
│   └── tables/
│
└── docs/
    └── thesis_summary.md
```

## Reproducing the Analysis

### Prerequisites

**R packages:**
```r
install.packages(c(
  "openair", "worldmet", "rmweather", "ranger",
  "augsynth", "dplyr", "tidyr", "purrr", "lubridate",
  "zoo", "imputeTS", "ggplot2", "readr", "stringr"
))
```

**Python packages** (for ERA5 download only):
```bash
pip install cdsapi xarray netcdf4 pandas numpy
```

You will also need a [CDS API key](https://cds.climate.copernicus.eu/api-how-to) configured in `~/.cdsapirc`.

### Run Order

1. `scripts/01_download_era5_weather.py` — Download ERA5 weather variables for all sites
2. `scripts/02_prepare_weather.R` — Merge NOAA surface weather with ERA5 variables
3. `scripts/03_import_and_normalize.R` — Import pollution data, impute gaps, run RF weather normalization
4. `scripts/04_build_panel.R` — Aggregate to weekly panel, assign treatment cohorts
5. `scripts/05_synthetic_control.R` — Select optimal nu, estimate staggered SCM
6. `scripts/06_robustness_checks.R` — Spatial/temporal placebo tests and LOOCV

## Technical Highlights

- **56 monitoring sites** across the UK (18 treated London sites + 38 control sites)
- **~61,000 hourly observations per site** (2018–2024)
- **Random Forest models** with OOB R² of 0.77–0.88, confirming strong weather-pollution relationships
- **Weather normalization reduces SD by 29–67%**, effectively removing meteorological noise
- **Optimal nu selection** using Ben-Michael et al. (2021) balance frontier heuristic
- **500-permutation spatial placebo tests** confirm treatment effects are not driven by chance

## References

- Ben-Michael, E., Feller, A., & Rothstein, J. (2021). Synthetic Controls with Staggered Adoption. *Journal of the Royal Statistical Society Series B*.
- Grange, S. K., & Carslaw, D. C. (2019). Using meteorological normalisation to detect interventions in air quality time series. *Science of the Total Environment*.
- Greater London Authority. ULEZ Expansion Reports.

## License

MIT License. See [LICENSE](LICENSE) for details.
