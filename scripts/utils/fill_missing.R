# ============================================================
# Missing Value Imputation Functions
# ============================================================
# Strategy:
# 1. Kalman Smoothing (StructTS) for short to medium gaps (<= 48 hours)
# 2. Gaps > 48 hours are left as NA to avoid artificial bias
# 3. Random Forest (rmweather) is robust to remaining NAs
# ============================================================

library(dplyr)
library(imputeTS)

#' Check missing value percentages for NO2 and NOx
#'
#' @param df Data frame with no2 and nox columns
#' @param site_name Character string for display
check_missing <- function(df, site_name) {
  missing_no2 <- sum(is.na(df$no2)) / nrow(df) * 100
  missing_nox <- sum(is.na(df$nox)) / nrow(df) * 100

  cat("Site:", site_name, "\n")
  cat("  NO2 missing:", round(missing_no2, 2), "%\n")
  cat("  NOx missing:", round(missing_nox, 2), "%\n")

  invisible(data.frame(
    site = site_name,
    no2_missing_pct = missing_no2,
    nox_missing_pct = missing_nox
  ))
}


#' Fill missing values using Kalman smoothing
#'
#' Applies Kalman smoothing (StructTS model) for gaps <= 48 hours.
#' Longer gaps are left as NA to avoid introducing bias.
#' Negative values are clipped to zero (physical constraint).
#'
#' @param df Data frame with date, no2, and nox columns
#' @param site_code Character string for the site ID
#' @param max_gap Integer, maximum gap length to impute (default: 48 hours)
#' @return Data frame with imputed values
fill_missing <- function(df, site_code, max_gap = 48) {

  cat("============================================================\n")
  cat(paste("Processing Site:", site_code, "\n"))
  cat("============================================================\n")

  # --- Initial diagnosis ---
  total_rows <- nrow(df)
  no2_missing_init <- sum(is.na(df$no2))
  nox_missing_init <- sum(is.na(df$nox))

  cat(paste("Total observations  :", total_rows, "\n"))
  cat(paste("Initial NO2 missing :", no2_missing_init,
            "(", round(no2_missing_init / total_rows * 100, 2), "%)\n"))
  cat(paste("Initial NOx missing :", nox_missing_init,
            "(", round(nox_missing_init / total_rows * 100, 2), "%)\n\n"))

  # Helper: longest continuous NA sequence
  get_max_gap <- function(vec) {
    if (sum(is.na(vec)) == 0) return(0)
    rl <- rle(is.na(vec))
    return(max(rl$lengths[rl$values]))
  }

  cat(paste("Max continuous gap (NO2):", get_max_gap(df$no2), "hours\n"))
  cat(paste("Max continuous gap (NOx):", get_max_gap(df$nox), "hours\n\n"))

  # --- Apply Kalman smoothing ---
  cat(paste0("Applying Kalman Smoothing for gaps <= ", max_gap, " hours...\n"))

  df_filled <- df %>%
    mutate(
      no2 = na_kalman(no2, model = "StructTS", maxgap = max_gap),
      nox = na_kalman(nox, model = "StructTS", maxgap = max_gap)
    ) %>%
    mutate(
      no2 = pmax(no2, 0),
      nox = pmax(nox, 0)
    )

  # --- Validation ---
  no2_missing_final <- sum(is.na(df_filled$no2))
  nox_missing_final <- sum(is.na(df_filled$nox))

  cat("\n--- Imputation Summary ---\n")
  cat(paste("NO2: Filled", no2_missing_init - no2_missing_final, "points.",
            "Remaining NA:", no2_missing_final, "\n"))
  cat(paste("NOx: Filled", nox_missing_init - nox_missing_final, "points.",
            "Remaining NA:", nox_missing_final, "\n"))

  if (no2_missing_final > (total_rows * 0.2)) {
    cat("\n[WARNING]: Significant gaps (>20%) remain at this site.\n")
  }

  cat("============================================================\n\n")

  return(df_filled)
}
