# ============================================================
# 06_robustness_checks.R
# Robustness checks for the staggered SCM results:
# 1. Spatial placebo tests (random treatment assignment)
# 2. Leave-one-out cross-validation
# ============================================================

library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(augsynth)

output_dir <- "outputs"

cat("========================================\n")
cat("ROBUSTNESS CHECKS\n")
cat("========================================\n\n")

# ============================================================
# Load data and previous results
# ============================================================

df_no2_model <- read_csv(file.path(output_dir, "panel_no2.csv"))
df_nox_model <- read_csv(file.path(output_dir, "panel_nox.csv"))
df_no2_model$week_date <- as.Date(df_no2_model$week_date)
df_nox_model$week_date <- as.Date(df_nox_model$week_date)

scm_results <- readRDS(file.path(output_dir, "scm_results.rds"))
optimal_nu_no2 <- scm_results$optimal_nu_no2
optimal_nu_nox <- scm_results$optimal_nu_nox
actual_att_no2 <- scm_results$actual_att_no2
actual_att_nox <- scm_results$actual_att_nox

# Treatment site list
treatment_sites <- c("EN1", "IS2", "CR5", "RI1", "HV3", "TH2", "GN3", "EA6",
                     "KC1", "BT8", "MY1", "WM6", "EI1", "LW4", "SK6", "EA8",
                     "WM5", "TH4")

# Cohort structure for placebo
cohort_structure <- data.frame(
  cohort = c("Central_2019", "Inner_2021", "Outer_2023"),
  treat_date = as.Date(c("2019-04-08", "2021-10-25", "2023-08-29")),
  n_sites = c(4, 7, 7),
  stringsAsFactors = FALSE
)


# ============================================================
# PART 1: SPATIAL PLACEBO TESTS
# ============================================================

run_spatial_placebo <- function(df_base, outcome_var, pollutant_name,
                                 actual_att, n_target = 500, max_attempts = 2000) {

  cat(paste0("\n[", pollutant_name, "] Running spatial placebo tests...\n"))
  cat(paste0("  Target: ", n_target, " successful permutations\n"))

  control_sites <- df_base %>%
    filter(is_treated == 0) %>%
    distinct(site_code) %>%
    pull(site_code)

  n_treat <- sum(cohort_structure$n_sites)
  placebo_results <- tibble()
  successful <- 0
  attempts <- 0

  while (successful < n_target && attempts < max_attempts) {
    attempts <- attempts + 1

    if (attempts %% 50 == 0) {
      cat(paste0("  Attempt ", attempts, " | Successful: ", successful, "\n"))
    }

    # Randomly assign control sites to fake cohorts
    shuffled <- sample(control_sites, min(n_treat, length(control_sites)))
    fake_cohort <- data.frame(
      site_code  = shuffled,
      treat_date = rep(cohort_structure$treat_date,
                       times = cohort_structure$n_sites)
    )

    df_spatial <- df_base %>%
      left_join(fake_cohort, by = "site_code") %>%
      mutate(
        placebo_treat = if_else(!is.na(treat_date) & week_date >= treat_date, 1, 0),
        treatment_date = if_else(!is.na(treat_date), treat_date, as.Date("2099-01-01"))
      ) %>%
      select(-treat_date)

    tryCatch({
      scm_placebo <- multisynth(
        as.formula(paste(outcome_var, "~ placebo_treat")),
        unit = site_code, time = week_date,
        data = df_spatial,
        n_lags = 52, n_leads = 52, nu = NULL, fixedeff = TRUE
      )

      att_placebo <- as.data.frame(summary(scm_placebo)$att) %>%
        filter(Level == "Average", Time >= 0, Time <= 52, !is.nan(Estimate)) %>%
        summarise(permutation = successful + 1,
                  mean_att = mean(Estimate))

      placebo_results <- bind_rows(placebo_results, att_placebo)
      successful <- successful + 1

      rm(scm_placebo)
      if (attempts %% 10 == 0) gc()

    }, error = function(e) NULL)
  }

  # Compute p-value
  cat(paste0("\n[", pollutant_name, "] Results:\n"))
  cat(paste0("  Total attempts:  ", attempts, "\n"))
  cat(paste0("  Successful:      ", successful, "\n"))
  cat(paste0("  Failure rate:    ", round(100 * (attempts - successful) / attempts, 1), "%\n"))

  if (nrow(placebo_results) > 0) {
    p_value <- mean(placebo_results$mean_att <= actual_att)

    cat(paste0("  Placebo ATT mean: ", round(mean(placebo_results$mean_att), 4), "\n"))
    cat(paste0("  Placebo ATT sd:   ", round(sd(placebo_results$mean_att), 4), "\n"))
    cat(paste0("  p-value:          ", round(p_value, 4), "\n"))

    # Plot
    p <- ggplot(placebo_results, aes(x = mean_att)) +
      geom_histogram(bins = 30, fill = "gray70", color = "black", alpha = 0.7) +
      geom_vline(xintercept = actual_att, color = "red", linewidth = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      annotate("text", x = actual_att - 0.005, y = Inf,
               label = paste0("Actual ATT = ", round(actual_att, 3)),
               color = "red", size = 5, hjust = 1, vjust = 2) +
      labs(
        title = paste0("Spatial Placebo Test: ", pollutant_name),
        subtitle = paste0(successful, " permutations | p = ", round(p_value, 3)),
        x = "Mean ATT (placebo)", y = "Count"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))

    ggsave(file.path(output_dir, "figures",
                     paste0("placebo_", tolower(pollutant_name), ".png")),
           p, width = 10, height = 6, dpi = 300)
    print(p)
  }

  return(placebo_results)
}

# Run placebo tests
# Use the base panel data (control sites only) for placebo permutations
df_no2_base <- df_no2_model %>% filter(is_treated == 0)
df_nox_base <- df_nox_model %>% filter(is_treated == 0)

placebo_no2 <- run_spatial_placebo(df_no2_base, "log_no2", "NO2", actual_att_no2)
placebo_nox <- run_spatial_placebo(df_nox_base, "log_nox", "NOx", actual_att_nox)


# ============================================================
# PART 2: LEAVE-ONE-OUT CROSS-VALIDATION
# ============================================================

run_loocv <- function(df_model, outcome_var, pollutant_name,
                       optimal_nu, actual_att) {

  cat(paste0("\n[", pollutant_name, "] Running leave-one-out CV...\n"))

  loocv_results <- tibble()

  for (site in treatment_sites) {
    cat(paste0("  Removing site: ", site, "\n"))

    df_loo <- df_model %>% filter(site_code != site)

    tryCatch({
      scm_loo <- multisynth(
        as.formula(paste(outcome_var, "~ is_treated")),
        unit = site_code, time = week_date,
        data = df_loo,
        n_lags = 52, n_leads = 52,
        nu = optimal_nu, fixedeff = TRUE
      )

      att_loo <- as.data.frame(summary(scm_loo)$att) %>%
        filter(Level == "Average", Time >= 0, Time <= 52, !is.nan(Estimate)) %>%
        summarise(
          removed_site = site,
          mean_att     = mean(Estimate),
          se           = sd(Estimate) / sqrt(n()),
          ci_lower     = mean_att - 1.96 * se,
          ci_upper     = mean_att + 1.96 * se
        )

      loocv_results <- bind_rows(loocv_results, att_loo)
    }, error = function(e) {
      cat(paste0("    Error removing ", site, ": ", e$message, "\n"))
    })
  }

  cat(paste0("\n[", pollutant_name, "] LOOCV Results:\n"))
  print(loocv_results %>% mutate(att_pct = round((exp(mean_att) - 1) * 100, 2)))

  # Plot
  p <- ggplot(loocv_results, aes(x = mean_att, y = reorder(removed_site, mean_att))) +
    geom_vline(xintercept = actual_att, linetype = "dashed", color = "red", size = 1) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
    geom_point(color = "steelblue", size = 3) +
    labs(
      title = paste0("Leave-One-Out Cross Validation: ", pollutant_name),
      subtitle = "Effect of removing each treated site",
      x = "Average ATT", y = "Removed Site"
    ) +
    theme_minimal()

  ggsave(file.path(output_dir, "figures",
                   paste0("loocv_", tolower(pollutant_name), ".png")),
         p, width = 10, height = 7, dpi = 300)
  print(p)

  return(loocv_results)
}

# Run LOOCV
loocv_no2 <- run_loocv(df_no2_model, "log_no2", "NO2", optimal_nu_no2, actual_att_no2)
loocv_nox <- run_loocv(df_nox_model, "log_nox", "NOx", optimal_nu_nox, actual_att_nox)

# ============================================================
# Save robustness check results
# ============================================================

write_csv(placebo_no2, file.path(output_dir, "placebo_no2.csv"))
write_csv(placebo_nox, file.path(output_dir, "placebo_nox.csv"))
write_csv(loocv_no2, file.path(output_dir, "loocv_no2.csv"))
write_csv(loocv_nox, file.path(output_dir, "loocv_nox.csv"))

cat("\n========================================\n")
cat("ROBUSTNESS CHECKS COMPLETE\n")
cat("========================================\n")
