# ============================================================
# Weather Normalization Analysis Function
# ============================================================
# Uses Random Forest (rmweather) to:
# 1. Model pollution as a function of weather + temporal variables
# 2. Generate weather-normalized concentrations by shuffling weather
# 3. Produce diagnostic plots and model performance metrics
#
# Supports KCL, AURN, NI, and WAQN data sources.
# ============================================================

library(rmweather)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(zoo)

#' Analyze a single monitoring site
#'
#' @param site_code Character, site identifier (e.g., "EN1")
#' @param site_name Character, human-readable site name
#' @param site_data Data frame with pre-processed pollution data (date, no2, nox)
#' @param weather_data Data frame with weather variables (optional)
#' @param weather_file Path to weather CSV file (optional)
#' @param source Character, data source: "kcl", "aurn", "ni", or "waqn"
#' @param output_dir Character, directory for output files (default: "outputs")
#' @return List with weekly data, model performance, and SD summaries
analyze_site <- function(site_code, site_name,
                         site_data,
                         weather_data = NULL,
                         weather_file = NULL,
                         source = "kcl",
                         output_dir = "outputs") {

  # ---- Academic plot theme ----
  theme_academic <- function() {
    theme_minimal() +
      theme(
        text             = element_text(size = 10),
        axis.text        = element_text(size = 9, color = "black"),
        axis.title       = element_text(size = 10),
        axis.line        = element_line(color = "black", size = 0.5),
        panel.grid.major = element_line(color = "grey90", size = 0.25),
        panel.grid.minor = element_blank(),
        legend.position  = "top",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 9),
        plot.title       = element_text(size = 11, hjust = 0.5),
        plot.subtitle    = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.margin      = margin(5, 5, 5, 5, "mm")
      )
  }

  col_obs    <- "#8B8B8B"
  col_norm   <- "#32CD32"
  col_policy <- "black"

  # ---- Helper: robust prediction data frame ----
  .predict_df <- function(model, df) {
    out <- tryCatch(rmw_predict(model = model, df = df), error = function(e) e)
    if (inherits(out, "error")) {
      if (!is.null(model$predictions) && length(model$predictions) >= nrow(df)) {
        tibble(date = df$date, predicted = as.numeric(model$predictions[seq_len(nrow(df))]))
      } else {
        stop("Could not get predictions.")
      }
    } else if (is.numeric(out)) {
      tibble(date = df$date, predicted = as.numeric(out))
    } else if (is.data.frame(out)) {
      if ("value_predict" %in% names(out) && "date" %in% names(out)) {
        tibble(date = out$date, predicted = out$value_predict)
      } else if ("value_predict" %in% names(out)) {
        tibble(date = df$date, predicted = out$value_predict)
      } else if ("predicted" %in% names(out)) {
        tibble(date = if ("date" %in% names(out)) out$date else df$date,
               predicted = out$predicted)
      } else {
        tibble(date = df$date, predicted = as.numeric(out[[1]]))
      }
    } else {
      tibble(date = df$date, predicted = as.numeric(out))
    }
  }

  .dedup_pred <- function(pred_df) {
    pred_df %>%
      group_by(date) %>%
      summarise(predicted = mean(predicted, na.rm = TRUE), .groups = "drop")
  }

  # ============================================================
  # Step 0: Load weather data
  # ============================================================
  if (!is.null(weather_data)) {
    weather_final <- weather_data
  } else if (!is.null(weather_file)) {
    weather_final <- read.csv(weather_file)
  } else {
    stop(paste("No weather data provided for", site_code))
  }

  # Parse date column
  if ("date" %in% names(weather_final)) {
    weather_final$date <- as.POSIXct(weather_final$date, tz = "UTC")
  } else if ("date_unix" %in% names(weather_final)) {
    weather_final$date <- as.POSIXct(weather_final$date_unix, origin = "1970-01-01", tz = "UTC")
  } else if ("datetime" %in% names(weather_final)) {
    weather_final$date <- as.POSIXct(weather_final$datetime, tz = "UTC")
  } else {
    stop("Weather data must contain a 'date', 'date_unix', or 'datetime' column.")
  }

  # Standardize column names
  weather_final <- weather_final %>%
    rename_with(~ case_when(
      . == "blh_m"        ~ "boundary_layer_height",
      . == "solar_wm2"    ~ "solar_rad",
      . == "cloud_cover"  ~ "cloud_cover_pct",
      . == "precip_total" ~ "precip",
      TRUE ~ .
    ))

  # ============================================================
  # Step 1: Validate site data
  # ============================================================
  cat("========================================\n")
  cat(paste("Analyzing:", site_code, "-", site_name, "\n"))
  cat(paste("Source:", toupper(source), "\n"))
  cat("========================================\n")

  site_data$date <- as.POSIXct(site_data$date, tz = "UTC")
  cat(paste("Input data rows:", nrow(site_data), "\n"))

  # ============================================================
  # Step 2: Merge pollution and weather data
  # ============================================================
  site_model_data <- site_data %>%
    filter(!is.na(no2) & !is.na(nox)) %>%
    left_join(weather_final, by = "date")

  required_weather <- c("ws", "wd", "air_temp", "RH", "atmos_pres")
  optional_weather <- c("solar_rad", "boundary_layer_height", "cloud_cover_pct", "precip")

  missing_optional <- setdiff(optional_weather, names(site_model_data))
  if (length(missing_optional) > 0) {
    for (var in missing_optional) site_model_data[[var]] <- NA
  }

  site_model_data <- site_model_data %>%
    filter(complete.cases(ws, wd, air_temp, RH, atmos_pres))

  # Policy dates for vertical lines on plots
  if (tolower(source) == "kcl") {
    policy_dates  <- as.POSIXct(c("2019-04-08", "2021-10-25", "2023-08-29"), tz = "UTC")
    policy_labels <- c("ULEZ Launch", "Expansion I", "Expansion II")
  } else if (tolower(source) == "waqn") {
    policy_dates  <- as.POSIXct(c("2018-06-18", "2020-03-01"), tz = "UTC")
    policy_labels <- c("50mph Zones", "Clean Air Plan")
  } else {
    policy_dates  <- as.POSIXct(character(0), tz = "UTC")
    policy_labels <- character(0)
  }

  # ============================================================
  # Step 3: NO2 - Random Forest + Weather Normalization
  # ============================================================
  cat("\n===== NO2 ANALYSIS =====\n")
  set.seed(123)
  site_prepared_no2 <- rmw_prepare_data(df = site_model_data, value = "no2", na.rm = TRUE)

  available_variables <- c(
    "ws", "wd", "air_temp", "RH", "atmos_pres", "solar_rad",
    "boundary_layer_height", "cloud_cover_pct", "precip",
    "date_unix", "day_julian", "week", "weekday", "hour", "month"
  )
  available_variables <- available_variables[available_variables %in% names(site_prepared_no2)]
  available_variables <- available_variables[
    sapply(site_prepared_no2[available_variables], function(x) !all(is.na(x)))
  ]

  cat("Predictor variables:", paste(available_variables, collapse = ", "), "\n")

  rf_model_no2 <- rmw_train_model(
    df = site_prepared_no2, variables = available_variables,
    n_trees = 500, mtry = 5, min_node_size = 5, verbose = FALSE
  )
  no2_oob_r2 <- rf_model_no2$r.squared
  cat("NO2 OOB R2:", round(no2_oob_r2, 3), "\n")

  weather_variables <- intersect(
    c("ws", "wd", "air_temp", "RH", "atmos_pres",
      "solar_rad", "boundary_layer_height", "cloud_cover_pct", "precip"),
    available_variables
  )

  site_normalised_no2 <- rmw_normalise(
    model = rf_model_no2, df = site_prepared_no2,
    variables = weather_variables, n_samples = 500, verbose = FALSE
  )

  site_comparison_no2 <- site_normalised_no2 %>%
    left_join(site_prepared_no2 %>% select(date, observed = value), by = "date") %>%
    rename(normalised = value_predict)

  no2_sd_observed   <- sd(site_comparison_no2$observed, na.rm = TRUE)
  no2_sd_normalised <- sd(site_comparison_no2$normalised, na.rm = TRUE)
  no2_sd_change     <- ((no2_sd_normalised - no2_sd_observed) / no2_sd_observed) * 100

  site_weekly_no2 <- site_comparison_no2 %>%
    mutate(week_date = floor_date(date, "week", week_start = 1)) %>%
    group_by(week_date) %>%
    summarise(observed_weekly   = mean(observed, na.rm = TRUE),
              normalised_weekly = mean(normalised, na.rm = TRUE),
              n_obs = n(), .groups = "drop") %>%
    filter(n_obs >= 100)

  pred_no2_hourly <- .predict_df(rf_model_no2, site_prepared_no2) %>% .dedup_pred()

  no2_monthly <- site_prepared_no2 %>%
    select(date, observed = value) %>%
    left_join(pred_no2_hourly, by = "date") %>%
    mutate(month_year = floor_date(date, "month")) %>%
    group_by(month_year) %>%
    summarise(observed_monthly  = mean(observed, na.rm = TRUE),
              predicted_monthly = mean(predicted, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    filter(n >= 500)

  no2_monthly_r2 <- suppressWarnings(
    cor(no2_monthly$observed_monthly, no2_monthly$predicted_monthly, use = "complete.obs")^2
  )

  # ============================================================
  # Step 4: NOx - Random Forest + Weather Normalization
  # ============================================================
  cat("\n===== NOx ANALYSIS =====\n")
  set.seed(123)
  site_prepared_nox <- rmw_prepare_data(df = site_model_data, value = "nox", na.rm = TRUE)

  rf_model_nox <- rmw_train_model(
    df = site_prepared_nox, variables = available_variables,
    n_trees = 500, mtry = 5, min_node_size = 5, verbose = FALSE
  )
  nox_oob_r2 <- rf_model_nox$r.squared
  cat("NOx OOB R2:", round(nox_oob_r2, 3), "\n")

  site_normalised_nox <- rmw_normalise(
    model = rf_model_nox, df = site_prepared_nox,
    variables = weather_variables, n_samples = 500, verbose = FALSE
  )

  site_comparison_nox <- site_normalised_nox %>%
    left_join(site_prepared_nox %>% select(date, observed = value), by = "date") %>%
    rename(normalised = value_predict)

  nox_sd_observed   <- sd(site_comparison_nox$observed, na.rm = TRUE)
  nox_sd_normalised <- sd(site_comparison_nox$normalised, na.rm = TRUE)
  nox_sd_change     <- ((nox_sd_normalised - nox_sd_observed) / nox_sd_observed) * 100

  site_weekly_nox <- site_comparison_nox %>%
    mutate(week_date = floor_date(date, "week", week_start = 1)) %>%
    group_by(week_date) %>%
    summarise(observed_weekly   = mean(observed, na.rm = TRUE),
              normalised_weekly = mean(normalised, na.rm = TRUE),
              n_obs = n(), .groups = "drop") %>%
    filter(n_obs >= 100)

  pred_nox_hourly <- .predict_df(rf_model_nox, site_prepared_nox) %>% .dedup_pred()

  nox_monthly <- site_prepared_nox %>%
    select(date, observed = value) %>%
    left_join(pred_nox_hourly, by = "date") %>%
    mutate(month_year = floor_date(date, "month")) %>%
    group_by(month_year) %>%
    summarise(observed_monthly  = mean(observed, na.rm = TRUE),
              predicted_monthly = mean(predicted, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    filter(n >= 500)

  nox_monthly_r2 <- suppressWarnings(
    cor(nox_monthly$observed_monthly, nox_monthly$predicted_monthly, use = "complete.obs")^2
  )

  # ============================================================
  # Step 5: Summary
  # ============================================================
  sd_summary <- data.frame(
    Site           = site_code,
    Site_Name      = site_name,
    Source         = toupper(source),
    Pollutant      = c("NO2", "NOx"),
    OOB_R2         = c(no2_oob_r2, nox_oob_r2),
    SD_Observed    = c(no2_sd_observed, nox_sd_observed),
    SD_Normalised  = c(no2_sd_normalised, nox_sd_normalised),
    SD_Change_Rate = c(no2_sd_change, nox_sd_change),
    stringsAsFactors = FALSE
  )

  cat("\n===== MODEL PERFORMANCE SUMMARY =====\n")
  cat("NO2 OOB R2:", round(no2_oob_r2, 3), "| Monthly R2:", round(no2_monthly_r2, 3), "\n")
  cat("NOx OOB R2:", round(nox_oob_r2, 3), "| Monthly R2:", round(nox_monthly_r2, 3), "\n")

  # ============================================================
  # Step 6: Save outputs
  # ============================================================
  prefix <- file.path(output_dir, paste0(site_code, "_", toupper(source)))

  write.csv(site_weekly_no2, paste0(prefix, "_weekly_NO2.csv"), row.names = FALSE)
  write.csv(site_weekly_nox, paste0(prefix, "_weekly_NOx.csv"), row.names = FALSE)
  write.csv(sd_summary, paste0(prefix, "_SD_summary.csv"), row.names = FALSE)

  # ============================================================
  # Step 7: Clean up memory
  # ============================================================
  final_output <- list(
    no2_weekly     = site_weekly_no2,
    nox_weekly     = site_weekly_nox,
    no2_monthly    = no2_monthly,
    nox_monthly    = nox_monthly,
    no2_oob_r2     = no2_oob_r2,
    nox_oob_r2     = nox_oob_r2,
    sd_summary     = sd_summary,
    no2_sd_change  = no2_sd_change,
    nox_sd_change  = nox_sd_change
  )

  rm(rf_model_no2, rf_model_nox, site_model_data,
     site_prepared_no2, site_prepared_nox,
     pred_no2_hourly, pred_nox_hourly,
     site_comparison_no2, site_comparison_nox,
     site_normalised_no2, site_normalised_nox)
  gc()

  cat("Done:", site_code, "(", toupper(source), ")\n")
  cat("========================================\n")

  return(final_output)
}
