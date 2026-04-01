# Data Directory

Raw data files are not tracked in this repository due to file size. Below is how to obtain each dataset.

## Pollution Data

London treatment sites are downloaded via `openair::importKCL()`. Control group sites use `openair::importAURN()`, `openair::importWAQN()`, and `openair::importNI()` depending on the network. See `config/sites_config.csv` for the full list of sites and their data sources.

## Weather Data

**NOAA surface observations (London Heathrow):**
Downloaded via `worldmet::importNOAA(code = "037720-99999", year = 2018:2024)`.

**ERA5 reanalysis (all sites):**
Downloaded via the CDS API using `scripts/01_download_era5_weather.py`. Variables include wind components, temperature, dewpoint, surface pressure, cloud cover, boundary layer height, precipitation, and surface solar radiation. Requires a CDS API key — see https://cds.climate.copernicus.eu/api-how-to.

## Site Metadata

`site_metadata/london_valid_sites.csv` contains the filtered list of London monitoring sites active from 2018–2024, exported from `openair::importMeta(source = "kcl")`.
