#### This script performs the analysis on the simulation output ####
library(tidyverse)
#### analysis helper function ####
get_performance <- function(
  results,
  aspect_var = NULL,
  cond_cols = NULL,
  method_regex,
  param_regex = "(phi|zeta)",
  grouping_regex = "(g1|g2)",
  stat_regex = "(est|se)"
) {
  custom_regex <- paste0(
    "^",
    method_regex,
    "_",
    param_regex,
    "\\d{1,2}_",
    grouping_regex,
    "_",
    stat_regex
  )
  # bring results in wide format: 1 row per parameter (e.g., phi11), group, and
  # statistic (est or se)
  results_long <- results |>
    pivot_longer(
      cols = matches(custom_regex),
      names_to = c("method", "parameter", "group", "stat"),
      names_sep = "_"
    ) |>
    # and then back into wide format with est and se in different columns
    pivot_wider(names_from = stat, values_from = value)

  # get population values for each parameter in each group
  pop_values <- results |>
    select(matches("^(phi|zeta)\\d{1,2}_g\\d+_pop$")) |>
    pivot_longer(everything(), names_to = "pop_name", values_to = "pop") |>
    separate(pop_name, into = c("parameter", "group", "suffix"), sep = "_") |>
    distinct(parameter, group, pop, .keep_all = FALSE)

  # add population values to previous object
  results_long <- results_long |>
    left_join(pop_values, by = c("parameter", "group"))

  # standard error depends on condition, so we compute the sd(est) for all
  # conditions, method, parameter, group
  se_values <- results_long |>
    group_by(across(all_of(cond_cols)), method, parameter, group) |>
    summarize(est_sd = sd(est), .groups = "drop")
  # add to results_long
  results_long <- results_long |>
    left_join(se_values, by = c(cond_cols, "method", "parameter", "group"))

  # compute difference between estimate and population value, CI,  whether
  # the  CI covers the population value, and the SE ratio per replication
  results_long <- results_long |>
    mutate(
      difference = est - pop,
      ci_lower = est - qnorm(.975) * se,
      ci_upper = est + qnorm(.975) * se,
      covered = (ci_lower <= pop) & (ci_upper >= pop),
      SER = se / est_sd
    )

  # compute outcomes (bias, RMSE, SE recovery, coverage), depending on aspect_var
  outcomes <- results_long |>
    group_by(across(all_of(aspect_var)), method, parameter, group) |>
    summarize(
      AB = mean(difference),
      RMSE = sqrt(mean((difference)^2)),
      mean_SER = mean(SER),
      COV = mean(covered),
      .groups = "drop"
    ) |>
    # pivot into long format (1 row per outcome)
    pivot_longer(AB:COV, names_to = "outcome") |>
    mutate(column_name = paste(parameter, group, sep = "_")) |>
    select(method, outcome, all_of(aspect_var), column_name, value) |>
    # pivot into wide format (1 column per parameter per group (e.g., phi11_g1))
    pivot_wider(names_from = column_name, values_from = value)

  return(outcomes)
}

#### Simulation 1 ####
# names of the condition columns:
cond_cols <- c("invariance_level", "pattern", "ss_n", "ss_t", "ss_ratio")

# read in the data
results <- read_csv(
  "Simulation 1/sim1.csv",
  col_types = cols(
    step1_single_warning_text = col_character(),
    step2_single_warning_text = col_character(),
    step3_single_warning_text = col_character(),
    step1_single_error_text = col_character(),
    step2_single_error_text = col_character(),
    step3_single_error_text = col_character(),
    step1_multi_warning_text = col_character(),
    step2_multi_warning_text = col_character(),
    step3_multi_warning_text = col_character(),
    step1_multi_error_text = col_character(),
    step2_multi_error_text = col_character(),
    step3_multi_error_text = col_character()
  )
) |>
  mutate(
    invariance_level = factor(
      invariance_level,
      c("full_strict", "full_scalar", "partial_scalar", "partial_metric")
    ),
    pattern = factor(pattern, levels = c("unidirectional", "mixed")),
    ss_ratio = factor(ss_ratio, levels = c("balanced", "unbalanced"))
  ) |>
  arrange(iteration) # sort by iteration

# find "problematic" iterations:
results <- results |>
  mutate(
    # outliers: SE greater than 1
    outlier_SE = if_any(
      c(
        single_phi11_g1_se,
        single_phi22_g1_se,
        single_phi12_g1_se,
        single_phi21_g1_se,
        single_phi11_g2_se,
        single_phi22_g2_se,
        single_phi12_g2_se,
        single_phi21_g2_se,
        multi_phi11_g1_se,
        multi_phi22_g1_se,
        multi_phi12_g1_se,
        multi_phi21_g1_se,
        multi_phi11_g2_se,
        multi_phi22_g2_se,
        multi_phi12_g2_se,
        multi_phi21_g2_se
      ),
      ~ .x > 1
    ),
    # missing SEs (i.e., those that are NA)
    NA_SE = if_any(
      c(
        single_phi11_g1_se,
        single_phi22_g1_se,
        single_phi12_g1_se,
        single_phi21_g1_se,
        single_phi11_g2_se,
        single_phi22_g2_se,
        single_phi12_g2_se,
        single_phi21_g2_se,
        multi_phi11_g1_se,
        multi_phi22_g1_se,
        multi_phi12_g1_se,
        multi_phi21_g1_se,
        multi_phi11_g2_se,
        multi_phi22_g2_se,
        multi_phi12_g2_se,
        multi_phi21_g2_se
      ),
      ~ is.na(.x)
    ),
    # variable that catches all warnings, errors, SE outliers, and missing SEs:
    reestimate = if_any(c(
      outlier_SE,
      NA_SE,
      step1_single_warning,
      step1_multi_warning,
      step3_single_error,
      step3_multi_error
    ))
  )

# investigate estimation problems:
results |>
  summarize(across(
    c(step1_single_warning:rerun_step1_multi),
    ~ sum(.x, na.rm = TRUE)
  )) |>
  print(width = Inf)
# step1 warnings: 8 single-group, 15 multi-group
# step3 error: 606 single-group, 531 multi-group
sum(results$step1_single_warning | results$step1_multi_warning)
# 23 data sets with step 1 warning
sum(results$step3_single_error | results$step3_multi_error)
# 900 data sets with step 3 error
sum(
  results$outlier_SE &
    !(results$step3_multi_error | results$step3_single_error),
  na.rm = TRUE
)
# 221 SEs greater than 1
sum(results$NA_SE & !(results$step3_multi_error | results$step3_single_error))
# 53 missing SEs

table(results$reestimate, results$invariance_level)
table(results$reestimate, results$pattern)
table(results$reestimate, results$ss_n)
table(results$reestimate, results$ss_t)
table(results$reestimate, results$ss_ratio)
# nearly all problems occured for small T (14 observations per person)
# also more problems when N is small (24 persons), and unbalanced sample sizes

# read in data from reestimation
results_reestimation <- read_csv(
  "Simulation 1/sim1_reestimation.csv",
  col_types = cols(
    step1_single_warning_text = col_character(),
    step2_single_warning_text = col_character(),
    step3_single_warning_text = col_character(),
    step1_single_error_text = col_character(),
    step2_single_error_text = col_character(),
    step3_single_error_text = col_character(),
    step1_multi_warning_text = col_character(),
    step2_multi_warning_text = col_character(),
    step3_multi_warning_text = col_character(),
    step1_multi_error_text = col_character(),
    step2_multi_error_text = col_character(),
    step3_multi_error_text = col_character()
  )
) |>
  mutate(
    invariance_level = factor(
      invariance_level,
      c("full_strict", "full_scalar", "partial_scalar", "partial_metric")
    ),
    pattern = factor(pattern, levels = c("unidirectional", "mixed")),
    ss_ratio = factor(ss_ratio, levels = c("balanced", "unbalanced"))
  ) |>
  arrange(iteration) # sort by iteration

# find "problematic" iterations:
results_reestimation <- results_reestimation |>
  mutate(
    # outliers: SE greater than 1
    outlier_SE = if_any(
      c(
        single_phi11_g1_se,
        single_phi22_g1_se,
        single_phi12_g1_se,
        single_phi21_g1_se,
        single_phi11_g2_se,
        single_phi22_g2_se,
        single_phi12_g2_se,
        single_phi21_g2_se,
        multi_phi11_g1_se,
        multi_phi22_g1_se,
        multi_phi12_g1_se,
        multi_phi21_g1_se,
        multi_phi11_g2_se,
        multi_phi22_g2_se,
        multi_phi12_g2_se,
        multi_phi21_g2_se
      ),
      ~ .x > 1
    ),
    # missing SEs (i.e., those that are NA)
    NA_SE = if_any(
      c(
        single_phi11_g1_se,
        single_phi22_g1_se,
        single_phi12_g1_se,
        single_phi21_g1_se,
        single_phi11_g2_se,
        single_phi22_g2_se,
        single_phi12_g2_se,
        single_phi21_g2_se,
        multi_phi11_g1_se,
        multi_phi22_g1_se,
        multi_phi12_g1_se,
        multi_phi21_g1_se,
        multi_phi11_g2_se,
        multi_phi22_g2_se,
        multi_phi12_g2_se,
        multi_phi21_g2_se
      ),
      ~ is.na(.x)
    ),
    # variable that catches all warnings, errors, SE outliers, and missing SEs:
    reestimate = if_any(c(
      outlier_SE,
      NA_SE,
      step1_single_warning,
      step1_multi_warning,
      step3_single_error,
      step3_multi_error
    ))
  )

# did the reestimation solve the problems?
results_reestimation |>
  summarize(across(
    c(step1_single_warning:rerun_step1_multi),
    ~ sum(.x, na.rm = TRUE)
  )) |>
  print(width = Inf)
# step1 warnings: 8 single-group, 15 multi-group
# step3 error: 449 (down from 606) single-group, 433 (down from 531) multi-group
sum(
  results_reestimation$step1_single_warning |
    results_reestimation$step1_multi_warning
)
# 23 data sets with step 1 warning
sum(
  results_reestimation$step3_single_error |
    results_reestimation$step3_multi_error
)
# 648 (down from 900) data sets with step 3 error
sum(
  results_reestimation$outlier_SE &
    !(results_reestimation$step3_multi_error |
      results_reestimation$step3_single_error),
  na.rm = TRUE
)
# 166 (down from 221) SEs greater than 1
sum(
  results_reestimation$NA_SE &
    !(results_reestimation$step3_multi_error |
      results_reestimation$step3_single_error)
)
# 33 (down from 53) missing SEs

# replace re-estimated data sets in results object:
results[
  results$iteration %in% results_reestimation$iteration,
] <- results_reestimation
rm(results_reestimation)

## compute convergence rates for each aspect
convergence_rates <- map(
  c(list(NULL), as.list(cond_cols)),
  function(aspect_var) {
    outcomes <- results |>
      group_by(across(all_of(aspect_var))) |>
      summarize(
        convergence_rate_single = mean(step3_single_error == FALSE),
        convergence_rate_multi = mean(step3_multi_error == FALSE),
        .groups = "drop"
      )
    # add "aspect" and "level" columns
    if (is.null(aspect_var)) {
      outcomes <- outcomes |>
        mutate(
          aspect = "overall",
          level = NA_character_,
          .before = "convergence_rate_single"
        )
    } else {
      outcomes <- outcomes |>
        rename(level = !!aspect_var) |>
        mutate(aspect = aspect_var, .before = level) |>
        mutate(level = as.character(level))
    }
  }
) |>
  list_rbind()

# write_csv(convergence_rates, file = "convergence_rates_study1.csv")

# remove data sets with SEs > 1, missing SEs, inadmissible solution in step 1,
# or failed convergence in step 3
results <- results |>
  mutate(
    remove = if_any(c(
      outlier_SE,
      NA_SE,
      step1_single_warning,
      step1_multi_warning,
      step3_single_error,
      step3_multi_error
    ))
  )
sum(results$remove, na.rm = TRUE)
# 860 data sets to remove

results <- results |>
  dplyr::filter(!remove)

# rename estimate columns (add _est suffix):
est_cols <- names(results)[
  str_detect(names(results), "^(single|multi)_(phi|zeta)\\d{1,2}_(g1|g2)$")
]
names(est_cols) <- paste0(est_cols, "_est")
results <- results |>
  rename(all_of(est_cols))


## recovery of MM parameters:
results |>
  summarize(
    mean_bias_lambda = mean(bias_lambda),
    mean_bias_theta = mean(bias_theta),
    mean_bias_tau = mean(bias_tau),
    mean_RMSE_lambda = mean(RMSE_lambda),
    mean_RMSE_theta = mean(RMSE_theta),
    mean_RMSE_tau = mean(RMSE_tau)
  )


## recovery of SM parameters
# compute marginal outcomes for levels of each aspect
# (averaged across all other aspects)
performance <- map(
  c(list(NULL), as.list(cond_cols)),
  function(aspect_var) {
    outcomes <- get_performance(
      results,
      aspect_var = aspect_var,
      method_regex = "(single|multi)",
      cond_cols = cond_cols
    )

    # add "aspect" and "level" columns
    if (is.null(aspect_var)) {
      outcomes <- outcomes |>
        mutate(aspect = "overall", level = NA_character_, .after = "outcome")
    } else {
      outcomes <- outcomes |>
        rename(level = !!aspect_var) |>
        mutate(aspect = aspect_var, .before = level) |>
        mutate(level = as.character(level))
    }
  }
) |>
  list_rbind()

# interaction between pattern and ss_ratio aspects in single-group
results |>
  get_performance(
    aspect_var = c("pattern", "ss_ratio"),
    method_regex = "single",
    cond_cols = cond_cols
  ) |>
  filter(outcome == "AB")

# create CSV files:
# write_csv(
#   performance |> dplyr::filter(outcome == "AB", method == "multi"),
#   file = "table_bias_multi.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "RMSE", method == "multi"),
#   file = "table_RMSE_multi.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "COV", method == "multi"),
#   file = "table_coverage_multi.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "mean_SER", method == "multi"),
#   file = "table_SER_multi.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "AB", method == "single"),
#   file = "table_bias_single.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "RMSE", method == "single"),
#   file = "table_RMSE_single.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "COV", method == "single"),
#   file = "table_coverage_single.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "mean_SER", method == "single"),
#   file = "table_SER_single.csv"
# )

# clean up:
rm(performance, results, cond_cols, est_cols)
#### Simulation 2 ####
# read in the data
results <- read_csv(
  "Simulation 2/sim2.csv",
  col_types = cols(
    step1_personspec_warning_text = col_character(),
    step2_personspec_warning_text = col_character(),
    step3_personspec_warning_text = col_character(),
    step1_personspec_error_text = col_character(),
    step2_personspec_error_text = col_character(),
    step3_personspec_error_text = col_character(),
    step1_groupspec_warning_text = col_character(),
    step2_groupspec_warning_text = col_character(),
    step3_groupspec_warning_text = col_character(),
    step1_groupspec_error_text = col_character(),
    step2_groupspec_error_text = col_character(),
    step3_groupspec_error_text = col_character()
  )
) |>
  mutate(
    heterogeneity = factor(heterogeneity, levels = c("small", "large")),
    nature = factor(nature, levels = c("continuous", "categorical"))
  ) |>
  arrange(iteration) # sort by iteration

## Warnings and errors:
results |>
  summarize(across(
    step1_personspec_warning:rerun_step1_groupspec,
    ~ sum(.x, na.rm = TRUE)
  )) |>
  print(width = Inf)
# warnings:
table(results$step1_personspec_warning_text)
# 13 convergence failures in step 1

sum(results$n_negvars_personspec > 0)
# 202 data sets where negative variances were estimated
table(results$n_negvars_personspec)
# 0    1    2    3
# 2198  179   22    1

# errors:
table(results$step3_groupspec_error_text)
# errors: 7 non-converged step3 models

# SEs greater than 1 or missing:
results <- results |>
  mutate(
    # outliers: SE greater than 1
    outlier_SE = if_any(
      c(
        personspec_phi11_g1_se,
        personspec_phi22_g1_se,
        personspec_phi12_g1_se,
        personspec_phi21_g1_se,
        personspec_phi11_g2_se,
        personspec_phi22_g2_se,
        personspec_phi12_g2_se,
        personspec_phi21_g2_se,
        groupspec_phi11_g1_se,
        groupspec_phi22_g1_se,
        groupspec_phi12_g1_se,
        groupspec_phi21_g1_se,
        groupspec_phi11_g2_se,
        groupspec_phi22_g2_se,
        groupspec_phi12_g2_se,
        groupspec_phi21_g2_se
      ),
      ~ .x > 1
    ),
    # missing SEs (i.e., those that are NA)
    NA_SE = if_any(
      c(
        personspec_phi11_g1_se,
        personspec_phi22_g1_se,
        personspec_phi12_g1_se,
        personspec_phi21_g1_se,
        personspec_phi11_g2_se,
        personspec_phi22_g2_se,
        personspec_phi12_g2_se,
        personspec_phi21_g2_se,
        groupspec_phi11_g1_se,
        groupspec_phi22_g1_se,
        groupspec_phi12_g1_se,
        groupspec_phi21_g1_se,
        groupspec_phi11_g2_se,
        groupspec_phi22_g2_se,
        groupspec_phi12_g2_se,
        groupspec_phi21_g2_se
      ),
      ~ is.na(.x)
    )
  )

sum(results$outlier_SE, na.rm = TRUE)
# 0 SEs > 1
sum(results$NA_SE)
# 7 missing SEs

results <- results |>
  mutate(
    remove = if_any(c(NA_SE, step1_personspec_warning, step3_groupspec_error))
  )
sum(results$remove, na.rm = TRUE)
# 20 data sets to remove

results <- results |>
  dplyr::filter(!remove)


# rename estimate columns (add _est suffix):
est_cols <- names(results)[
  str_detect(
    names(results),
    "^(personspec|groupspec)_(phi|zeta)\\d{1,2}_(g1|g2)$"
  )
]

names(est_cols) <- paste0(est_cols, "_est")

results <- results |>
  rename(all_of(est_cols))

## recovery of MM parameters (for person-specific models):
results |>
  summarize(
    mean_bias_lambda3 = mean(bias_lambda3),
    mean_bias_theta1 = mean(bias_theta1),
    mean_bias_theta3 = mean(bias_theta3),
    mean_RMSE_lambda3 = mean(RMSE_lambda3),
    mean_RMSE_theta1 = mean(RMSE_theta1),
    mean_RMSE_theta3 = mean(RMSE_theta3),
  ) |>
  print(width = Inf)


## recovery of SM parameters:
cond_cols <- c("ss_t", "heterogeneity", "nature")
performance <- map(
  c(list(NULL), as.list(cond_cols)),
  ~ get_performance(
    results,
    aspect_var = .x,
    method_regex = "(personspec|groupspec)",
    cond_cols = cond_cols
  )
) |>
  list_rbind()

# write_csv(
#   performance |> dplyr::filter(outcome == "AB", method == "personspec"),
#   file = "table_bias_personspec.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "RMSE", method == "personspec"),
#   file = "table_RMSE_personspec.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "COV", method == "personspec"),
#   file = "table_coverage_personspec.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "mean_SER", method == "personspec"),
#   file = "table_SER_personspec.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "AB", method == "groupspec"),
#   file = "table_bias_groupspec.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "RMSE", method == "groupspec"),
#   file = "table_RMSE_groupspec.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "COV", method == "groupspec"),
#   file = "table_coverage_groupspec.csv"
# )
# write_csv(
#   performance |> dplyr::filter(outcome == "mean_SER", method == "groupspec"),
#   file = "table_SER_groupspec.csv"
# )
