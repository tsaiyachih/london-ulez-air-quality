# ============================================================
# 04_build_panel.R
# Consolidate weather-normalized weekly data into a panel
# dataset for synthetic control estimation.
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(lubridate)

# ============================================================
# Configuration
# ============================================================

sites <- read_csv("config/sites_config.csv")
output_dir <- "outputs"

# ULEZ policy dates
date_central <- as.Date("2019-04-08")
date_inner   <- as.Date("2021-10-25")
date_outer   <- as.Date("2023-08-29")

# Treatment cohort site codes
sites_central <- c("MY1", "WM6", "SK6", "WM5")
sites_inner   <- c("IS2", "RI1", "TH2", "KC1", "BT8", "LW4", "TH4")
sites_outer   <- c("EN1", "CR5", "HV3", "GN3", "EA6", "EI1", "EA8")

cat("========================================\n")
cat("BUILDING WEEKLY PANEL DATA\n")
cat("========================================\n\n")

# ============================================================
# Load weekly CSV files produced by 03_import_and_normalize.R
# ============================================================

load_weekly_files <- function(output_dir, pollutant = "NO2") {
  pattern <- paste0("_weekly_", pollutant, "\\.csv$")
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop(paste("No weekly", pollutant, "files found in", output_dir))
  }

  map_df(files, function(f) {
    site_code <- str_split(basename(f), "_")[[1]][1]
    read_csv(f, show_col_types = FALSE) %>%
      mutate(site_code = site_code,
             week_date = as.Date(week_date))
  })
}

df_no2_raw <- load_weekly_files(output_dir, "NO2")
df_nox_raw <- load_weekly_files(output_dir, "NOx")

cat("Loaded NO2 weekly data:", nrow(df_no2_raw), "rows,",
    n_distinct(df_no2_raw$site_code), "sites\n")
cat("Loaded NOx weekly data:", nrow(df_nox_raw), "rows,",
    n_distinct(df_nox_raw$site_code), "sites\n")

# ============================================================
# Assign treatment information
# ============================================================

assign_treatment_info <- function(df) {
  df %>%
    left_join(sites %>% select(site_code, group, cohort, treatment_date),
              by = "site_code") %>%
    mutate(
      treatment_date = as.Date(treatment_date),
      treatment_date = if_else(group == "Control", as.Date("2099-01-01"), treatment_date),
      is_treated = if_else(week_date >= treatment_date, 1L, 0L)
    )
}

df_no2 <- assign_treatment_info(df_no2_raw)
df_nox <- assign_treatment_info(df_nox_raw)

cat("\nCohort distribution:\n")
print(table(df_no2 %>% distinct(site_code, cohort) %>% pull(cohort)))

# ============================================================
# Define donor pool (pure control sites)
# ============================================================

control_sites <- sites %>% filter(group == "Control") %>% pull(site_code)

df_no2_model <- df_no2 %>%
  filter(site_code %in% control_sites |
         (group == "Treatment" & week_date >= (treatment_date - 365)))

df_nox_model <- df_nox %>%
  filter(site_code %in% control_sites |
         (group == "Treatment" & week_date >= (treatment_date - 365)))

# Log transformation
df_no2_model <- df_no2_model %>% mutate(log_no2 = log(normalised_weekly + 1))
df_nox_model <- df_nox_model %>% mutate(log_nox = log(normalised_weekly + 1))

cat("\nModel data ready:\n")
cat("  NO2:", nrow(df_no2_model), "obs,", n_distinct(df_no2_model$site_code), "sites\n")
cat("  NOx:", nrow(df_nox_model), "obs,", n_distinct(df_nox_model$site_code), "sites\n")

# ============================================================
# Save panel datasets
# ============================================================

write_csv(df_no2_model, file.path(output_dir, "panel_no2.csv"))
write_csv(df_nox_model, file.path(output_dir, "panel_nox.csv"))

cat("\nSaved:\n")
cat("  outputs/panel_no2.csv\n")
cat("  outputs/panel_nox.csv\n")

cat("\n========================================\n")
cat("PANEL DATA COMPLETE\n")
cat("========================================\n")
