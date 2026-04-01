# ============================================================
# 02_prepare_weather.R
# Merge NOAA surface weather with ERA5 reanalysis data
# and handle missing values via linear interpolation.
# ============================================================

library(dplyr)
library(readr)
library(lubridate)
library(zoo)
library(worldmet)

# ============================================================
# London: NOAA Heathrow + ERA5 boundary layer height & solar
# ============================================================

cat("========================================\n")
cat("Preparing London weather data\n")
cat("========================================\n")

# Import NOAA surface observations from Heathrow Airport
heathrow_noaa <- importNOAA(
  code = "037720-99999",
  year = 2018:2024,
  hourly = TRUE,
  n.cores = 1
)

heathrow_noaa <- heathrow_noaa %>%
  mutate(
    date_unix    = as.numeric(date),
    year         = year(date),
    month        = month(date),
    week         = week(date),
    weekday      = wday(date),
    hour         = hour(date),
    day_julian   = yday(date),
    cloud_cover  = cl,
    precip_total = coalesce(precip, precip_6, precip_12, 0)
  ) %>%
  select(date, ws, wd, air_temp, RH, atmos_pres,
         cloud_cover, precip_total,
         date_unix, year, month, week, weekday, hour, day_julian)

# Import ERA5 variables (BLH and solar radiation for London)
# This file is produced by: python 01_download_era5_weather.py --city london ...
heathrow_era5 <- read_csv("data/raw/heathrow_blh_solar_processed.csv") %>%
  rename(date = datetime,
         blh_m = heathrow_blh_m,
         solar_wm2 = heathrow_solar_wm2) %>%
  mutate(date = as.POSIXct(date, tz = "UTC"))

# Merge NOAA + ERA5
london_weather_final <- heathrow_noaa %>%
  left_join(heathrow_era5, by = "date") %>%
  select(
    date, date_unix, year, month, week, weekday, hour, day_julian,
    ws, wd, air_temp, RH, atmos_pres, cloud_cover, precip_total,
    blh_m, solar_wm2
  )

# Check missing values
missing_check <- data.frame(
  missing_count = colSums(is.na(london_weather_final)),
  missing_pct   = colSums(is.na(london_weather_final)) / nrow(london_weather_final) * 100
)
cat("\nMissing values before interpolation:\n")
print(missing_check)

# Linear interpolation for small gaps (all < 2%)
cols_to_interpolate <- c("ws", "wd", "air_temp", "RH", "atmos_pres", "cloud_cover")
london_weather_final <- london_weather_final %>%
  mutate(across(all_of(cols_to_interpolate), ~ na.approx(., na.rm = FALSE)))

cat("\nMissing values after interpolation:\n")
missing_check_after <- data.frame(
  missing_count = colSums(is.na(london_weather_final)),
  missing_pct   = colSums(is.na(london_weather_final)) / nrow(london_weather_final) * 100
)
print(missing_check_after)

write_csv(london_weather_final, "data/processed/london_weather_final.csv")
cat("\nSaved: data/processed/london_weather_final.csv\n")


# ============================================================
# Control group cities: ERA5 only
# ============================================================
# For control group cities, all weather variables come from ERA5.
# Each city's weather file is produced by 01_download_era5_weather.py
# and already contains all required variables.
#
# Example:
#   python 01_download_era5_weather.py --city manchester --lat 53.3537 --lon -2.2749
#   python 01_download_era5_weather.py --city leeds --lat 53.8655 --lon -1.6609
#
# The output CSV files (e.g., manchester_weather_final.csv) can be
# used directly in 03_import_and_normalize.R via the weather_file parameter.

cat("\n========================================\n")
cat("Weather data preparation complete.\n")
cat("========================================\n")
