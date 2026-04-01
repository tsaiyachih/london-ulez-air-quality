# ============================================================
# 05_synthetic_control.R
# Partially Pooled Synthetic Control Method
# with Staggered Adoption (Ben-Michael et al., 2021)
#
# 1. Optimal nu selection via balance frontier
# 2. Main model estimation for NO2 and NOx
# 3. ATT extraction and visualization
# ============================================================

library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(augsynth)

output_dir <- "outputs"

cat("========================================\n")
cat("STAGGERED SYNTHETIC CONTROL ESTIMATION\n")
cat("========================================\n\n")

# ============================================================
# Load panel data
# ============================================================

df_no2_model <- read_csv(file.path(output_dir, "panel_no2.csv"))
df_nox_model <- read_csv(file.path(output_dir, "panel_nox.csv"))

df_no2_model$week_date <- as.Date(df_no2_model$week_date)
df_nox_model$week_date <- as.Date(df_nox_model$week_date)

# ============================================================
# Universal Nu Selection Function
# ============================================================

select_optimal_nu <- function(data, outcome_var, pollutant_name,
                               nu_grid = seq(0, 1, by = 0.05)) {

  cat("========================================\n")
  cat(paste0("OPTIMAL NU SELECTION - ", pollutant_name, "\n"))
  cat("Based on Ben-Michael et al. (2021)\n")
  cat("========================================\n\n")

  # Step 1: Compute balance frontier
  cat("[Step 1] Computing balance possibility frontier...\n")

  balance_results <- map_df(nu_grid, function(nu_val) {
    temp_mod <- multisynth(
      as.formula(paste(outcome_var, "~ is_treated")),
      unit = site_code, time = week_date,
      data = data, n_lags = 52, n_leads = 52, nu = nu_val
    )
    tibble(nu = nu_val, qpool = temp_mod$global_l2, qsep = temp_mod$ind_l2)
  })

  # Step 2: Heuristic nu
  L <- 52
  qpool_0 <- balance_results$qpool[balance_results$nu == 0]
  qsep_0  <- balance_results$qsep[balance_results$nu == 0]
  nu_heuristic <- min(sqrt(L) * qpool_0 / qsep_0, 1)

  cat(paste("  Heuristic nu:", round(nu_heuristic, 2), "\n"))

  # Step 3: Elbow point detection
  balance_results <- balance_results %>%
    mutate(qpool_improvement = lag(qpool) - qpool,
           qsep_cost = qsep - lag(qsep))

  initial_improvement <- balance_results$qpool_improvement[2]
  threshold <- 0.05 * initial_improvement

  elbow_candidates <- balance_results %>%
    filter(qpool_improvement < threshold, nu > 0) %>%
    arrange(nu)

  nu_elbow <- if (nrow(elbow_candidates) > 0) elbow_candidates$nu[1] else 0.5
  cat(paste("  Elbow nu:", round(nu_elbow, 2), "\n"))

  # Step 4: Test candidates with fixed effects
  cat("\n[Step 4] Testing candidates with fixed effects...\n")

  nu_candidates <- unique(c(0.5, nu_elbow, nu_heuristic, 0.6, 0.7, 0.8))
  nu_candidates <- sort(nu_candidates[nu_candidates >= 0 & nu_candidates <= 1])

  performance <- map_df(nu_candidates, function(nu_val) {
    test_mod <- multisynth(
      as.formula(paste(outcome_var, "~ is_treated")),
      unit = site_code, time = week_date,
      data = data, n_lags = 52, n_leads = 52,
      nu = nu_val, fixedeff = TRUE
    )
    summ <- summary(test_mod)
    tibble(
      nu = nu_val,
      scaled_qpool = summ$scaled_global_l2,
      scaled_qsep  = summ$scaled_ind_l2,
      avg_att = summ$att %>%
        filter(Level == "Average", Time >= 0) %>%
        summarise(mean(Estimate)) %>% pull()
    )
  })

  print(performance)

  # Step 5: Select optimal nu
  optimal <- performance %>%
    filter(scaled_qpool < 0.15) %>%
    arrange(scaled_qsep) %>%
    slice(1)

  if (nrow(optimal) == 0) {
    optimal <- performance %>%
      mutate(total_imbalance = scaled_qpool + scaled_qsep) %>%
      arrange(total_imbalance) %>%
      slice(1)
  }

  nu_optimal <- optimal$nu

  cat(paste("\nSelected nu:", round(nu_optimal, 2), "\n"))
  cat(paste("  Scaled Global L2:", round(optimal$scaled_qpool, 3), "\n"))
  cat(paste("  Scaled Individual L2:", round(optimal$scaled_qsep, 3), "\n"))

  # Step 6: Plot frontier
  p <- ggplot(balance_results, aes(x = qsep, y = qpool)) +
    geom_line(size = 1.2, color = "orange") +
    geom_point(size = 2) +
    geom_text(aes(label = round(nu, 2)), vjust = -0.5, size = 2) +
    labs(
      title = paste0(pollutant_name, ": Balance Possibility Frontier"),
      subtitle = paste0("Recommended nu = ", round(nu_optimal, 2)),
      x = "Unit-specific Imbalance (qsep)",
      y = "Pooled Imbalance (qpool)"
    ) +
    theme_minimal()

  ggsave(file.path(output_dir, "figures",
                   paste0("nu_frontier_", tolower(pollutant_name), ".png")),
         p, width = 8, height = 6, dpi = 300)
  print(p)

  list(nu_optimal = nu_optimal, nu_heuristic = nu_heuristic,
       nu_elbow = nu_elbow, balance_frontier = balance_results,
       performance_table = performance)
}


# ============================================================
# NO2 Analysis
# ============================================================

cat("\n\n========================================\n")
cat("NO2 ANALYSIS\n")
cat("========================================\n\n")

# Optimal nu selection
nu_result_no2 <- select_optimal_nu(df_no2_model, "log_no2", "NO2")
optimal_nu_no2 <- nu_result_no2$nu_optimal

# Main model
scm_no2 <- multisynth(
  log_no2 ~ is_treated,
  unit = site_code, time = week_date,
  data = df_no2_model,
  n_lags = 52, n_leads = 52,
  nu = optimal_nu_no2, fixedeff = TRUE
)

cat("\nNO2 Main Model Summary:\n")
print(summary(scm_no2))

# Extract ATT
att_no2 <- as.data.frame(summary(scm_no2)$att) %>%
  filter(Time >= -52, Time <= 52)

actual_att_no2 <- att_no2 %>%
  filter(Level == "Average", Time >= 0, Time <= 52, !is.nan(Estimate)) %>%
  summarise(mean_att = mean(Estimate)) %>%
  pull(mean_att)

cat(paste("NO2 Average ATT:", round(actual_att_no2, 4), "\n"))
cat(paste("NO2 Percentage change:", round((exp(actual_att_no2) - 1) * 100, 1), "%\n"))


# ============================================================
# NOx Analysis
# ============================================================

cat("\n\n========================================\n")
cat("NOx ANALYSIS\n")
cat("========================================\n\n")

# Optimal nu selection
nu_result_nox <- select_optimal_nu(df_nox_model, "log_nox", "NOx")
optimal_nu_nox <- nu_result_nox$nu_optimal

# Main model
scm_nox <- multisynth(
  log_nox ~ is_treated,
  unit = site_code, time = week_date,
  data = df_nox_model,
  n_lags = 52, n_leads = 52,
  nu = optimal_nu_nox, fixedeff = TRUE
)

cat("\nNOx Main Model Summary:\n")
print(summary(scm_nox))

# Extract ATT
att_nox <- as.data.frame(summary(scm_nox)$att) %>%
  filter(Time >= -52, Time <= 52)

actual_att_nox <- att_nox %>%
  filter(Level == "Average", Time >= 0, Time <= 52, !is.nan(Estimate)) %>%
  summarise(mean_att = mean(Estimate)) %>%
  pull(mean_att)

cat(paste("NOx Average ATT:", round(actual_att_nox, 4), "\n"))
cat(paste("NOx Percentage change:", round((exp(actual_att_nox) - 1) * 100, 1), "%\n"))


# ============================================================
# Save results
# ============================================================

write_csv(att_no2, file.path(output_dir, "att_no2.csv"))
write_csv(att_nox, file.path(output_dir, "att_nox.csv"))

saveRDS(list(scm_no2 = scm_no2, scm_nox = scm_nox,
             optimal_nu_no2 = optimal_nu_no2, optimal_nu_nox = optimal_nu_nox,
             actual_att_no2 = actual_att_no2, actual_att_nox = actual_att_nox),
        file.path(output_dir, "scm_results.rds"))

cat("\n========================================\n")
cat("SYNTHETIC CONTROL ESTIMATION COMPLETE\n")
cat("========================================\n")
