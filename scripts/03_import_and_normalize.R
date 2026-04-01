# ============================================================
# 03_import_and_normalize.R
# Import pollution data from KCL/AURN/WAQN/NI APIs,
# impute missing values, and run RF weather normalization.
#
# Driven by config/sites_config.csv — no need to copy-paste
# per-site code blocks.
# ============================================================

library(dplyr)
library(readr)
library(lubridate)
library(zoo)
library(openair)

# Load utility functions
source("scripts/utils/fill_missing.R")
source("scripts/utils/analyze_site.R")

# ============================================================
# Configuration
# ============================================================

sites <- read_csv("config/sites_config.csv")
weather_dir <- "data/processed"   # Directory containing *_weather_final.csv files
output_dir  <- "outputs"

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("POLLUTION DATA IMPORT & NORMALIZATION\n")
cat(paste("Total sites:", nrow(sites), "\n"))
cat("========================================\n\n")

# ============================================================
# Helper: Import pollution data based on source
# ============================================================

import_pollution <- function(site_code, source, years = 2018:2024) {
  src <- tolower(source)

  if (src == "kcl") {
    df <- importKCL(site = site_code, year = years)
  } else if (src == "aurn") {
    df <- importAURN(site = site_code, year = years) %>%
      select(date, no2, nox)
  } else if (src == "waqn") {
    df <- importWAQN(site = site_code, year = years) %>%
      select(date, no2, nox)
  } else if (src == "ni") {
    df <- importNI(site = site_code, year = years) %>%
      select(date, no2, nox)
  } else {
    stop(paste("Unknown source:", source))
  }

  df %>% mutate(date = with_tz(date, tzone = "UTC"))
}

# ============================================================
# Main loop: process each site
# ============================================================

results <- list()

for (i in seq_len(nrow(sites))) {
  site <- sites[i, ]

  cat("\n========================================\n")
  cat(paste0("[", i, "/", nrow(sites), "] ", site$site_code, " - ", site$site_name, "\n"))
  cat("========================================\n")

  # Step 1: Import pollution data
  tryCatch({
    raw <- import_pollution(site$site_code, site$source)

    # Step 2: Check and fill missing values
    check_missing(raw, site$site_code)
    filled <- fill_missing(raw, site$site_code)

    # Step 3: Load weather data
    weather_file <- file.path(weather_dir, paste0(site$weather_source, "_weather_final.csv"))
    if (!file.exists(weather_file)) {
      cat("  WARNING: Weather file not found:", weather_file, "\n")
      cat("  Skipping this site.\n")
      next
    }
    weather <- read.csv(weather_file)

    # Step 4: Run weather normalization
    results[[site$site_code]] <- analyze_site(
      site_code    = site$site_code,
      site_name    = site$site_name,
      site_data    = filled,
      weather_data = weather,
      source       = site$source,
      output_dir   = output_dir
    )

  }, error = function(e) {
    cat("  ERROR processing", site$site_code, ":", e$message, "\n")
  })
}

# ============================================================
# Summary
# ============================================================

cat("\n========================================\n")
cat("NORMALIZATION COMPLETE\n")
cat(paste("Successfully processed:", length(results), "/", nrow(sites), "sites\n"))
cat("========================================\n")

# Save results list
saveRDS(results, file.path(output_dir, "all_site_results.rds"))
cat("Saved: outputs/all_site_results.rds\n")
