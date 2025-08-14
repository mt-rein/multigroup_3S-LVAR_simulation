#### This script performs the analysis on the simulation output ####
library(tidyverse)
#### Simulation 1 ####
# read in the data
results <- read_csv("Simulation 1/sim1.csv",
                    col_types = cols(step1_single_warning_text = col_character(),
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
                                     step3_multi_error_text = col_character())) |> 
  mutate(invariance_level = factor(invariance_level, c("full_strict", "full_scalar", "partial_scalar", "partial_metric")),
         pattern = factor(pattern, levels = c("unidirectional", "mixed")),
         ss_ratio = factor(ss_ratio, levels = c("balanced", "unbalanced"))
  ) |> 
  arrange(iteration)                                                            # sort by iteration

# find "problematic" iterations:
results <- results |> 
  mutate(
    # outliers: SE greater than 1
    outlier_SE = if_any(c(
      single_phi11_g1_se, single_phi22_g1_se, single_phi12_g1_se, single_phi21_g1_se,
      single_phi11_g2_se, single_phi22_g2_se, single_phi12_g2_se, single_phi21_g2_se,
      multi_phi11_g1_se,  multi_phi22_g1_se,  multi_phi12_g1_se,  multi_phi21_g1_se,
      multi_phi11_g2_se,  multi_phi22_g2_se,  multi_phi12_g2_se,  multi_phi21_g2_se
    ), ~ .x > 1),
    # missing SEs (i.e., those that are NA)
    NA_SE = if_any(c(
      single_phi11_g1_se, single_phi22_g1_se, single_phi12_g1_se, single_phi21_g1_se,
      single_phi11_g2_se, single_phi22_g2_se, single_phi12_g2_se, single_phi21_g2_se,
      multi_phi11_g1_se,  multi_phi22_g1_se,  multi_phi12_g1_se,  multi_phi21_g1_se,
      multi_phi11_g2_se,  multi_phi22_g2_se,  multi_phi12_g2_se,  multi_phi21_g2_se
    ), ~ is.na(.x)),
    # variable that catches all warnings, errors, SE outliers, and missing SEs:
    reestimate = if_any(c(outlier_SE, NA_SE, step1_single_warning, step1_multi_warning, step3_single_error, step3_multi_error)))

# investigate estimation problems:
results |> 
  summarize(across(c(step1_single_warning:rerun_step1_multi), ~sum(.x, na.rm = TRUE))) |> 
  print(width = Inf)
# step1 warnings: 8 single-group, 15 multi-group
# step3 error: 606 single-group, 531 multi-group
sum(results$step1_single_warning | results$step1_multi_warning)
# 23 data sets with step 1 warning
sum(results$step3_single_error | results$step3_multi_error)
# 900 data sets with step 3 error
sum(results$outlier_SE & !(results$step3_multi_error | results$step3_single_error), na.rm = TRUE)
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
results_reestimation <- read_csv("Simulation 1/sim1_reestimation.csv",
                    col_types = cols(step1_single_warning_text = col_character(),
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
                                     step3_multi_error_text = col_character())) |> 
  mutate(invariance_level = factor(invariance_level, c("full_strict", "full_scalar", "partial_scalar", "partial_metric")),
         pattern = factor(pattern, levels = c("unidirectional", "mixed")),
         ss_ratio = factor(ss_ratio, levels = c("balanced", "unbalanced"))
  ) |> 
  arrange(iteration)                                                            # sort by iteration

# find "problematic" iterations:
results_reestimation <- results_reestimation |> 
  mutate(
    # outliers: SE greater than 1
    outlier_SE = if_any(c(
      single_phi11_g1_se, single_phi22_g1_se, single_phi12_g1_se, single_phi21_g1_se,
      single_phi11_g2_se, single_phi22_g2_se, single_phi12_g2_se, single_phi21_g2_se,
      multi_phi11_g1_se,  multi_phi22_g1_se,  multi_phi12_g1_se,  multi_phi21_g1_se,
      multi_phi11_g2_se,  multi_phi22_g2_se,  multi_phi12_g2_se,  multi_phi21_g2_se
    ), ~ .x > 1),
    # missing SEs (i.e., those that are NA)
    NA_SE = if_any(c(
      single_phi11_g1_se, single_phi22_g1_se, single_phi12_g1_se, single_phi21_g1_se,
      single_phi11_g2_se, single_phi22_g2_se, single_phi12_g2_se, single_phi21_g2_se,
      multi_phi11_g1_se,  multi_phi22_g1_se,  multi_phi12_g1_se,  multi_phi21_g1_se,
      multi_phi11_g2_se,  multi_phi22_g2_se,  multi_phi12_g2_se,  multi_phi21_g2_se
    ), ~ is.na(.x)),
    # variable that catches all warnings, errors, SE outliers, and missing SEs:
    reestimate = if_any(c(outlier_SE, NA_SE, step1_single_warning, step1_multi_warning, step3_single_error, step3_multi_error)))

# did the reestimation solve the problems?
results_reestimation |> 
  summarize(across(c(step1_single_warning:rerun_step1_multi), ~sum(.x, na.rm = TRUE))) |> 
  print(width = Inf)
# step1 warnings: 8 single-group, 15 multi-group
# step3 error: 449 (down from 606) single-group, 433 (down from 531) multi-group
sum(results_reestimation$step1_single_warning | results_reestimation$step1_multi_warning)
# 23 data sets with step 1 warning
sum(results_reestimation$step3_single_error | results_reestimation$step3_multi_error)
# 648 (down from 900) data sets with step 3 error
sum(results_reestimation$outlier_SE & !(results_reestimation$step3_multi_error | results_reestimation$step3_single_error), na.rm = TRUE)
# 166 (down from 221) SEs greater than 1
sum(results_reestimation$NA_SE & !(results_reestimation$step3_multi_error | results_reestimation$step3_single_error))
# 33 (down from 53) missing SEs

# replace re-estimated data sets in results object:
results[results$iteration %in% results_reestimation$iteration,] <- results_reestimation
rm(results_reestimation)

# remove data sets with SEs > 1, missing SEs, inadmissible solution in step 1, or failed convergence in step 3
results <- results |> 
  mutate(remove = if_any(c(outlier_SE, NA_SE, step1_single_warning, step1_multi_warning, step3_single_error, step3_multi_error)))
sum(results$remove, na.rm = TRUE)
# 860 data sets to remove

results <- results |> 
  dplyr::filter(!remove)

## recovery of MM parameters:
results |> group_by(ss_n, ss_t) |> 
  summarize(mean_bias_lambda = mean(bias_lambda),
            mean_bias_theta = mean(bias_theta),
            mean_bias_tau = mean(bias_tau), 
            mean_RMSE_lambda = mean(RMSE_lambda),
            mean_RMSE_theta = mean(RMSE_theta),
            mean_RMSE_tau = mean(RMSE_tau))

## recovery of SM parameters:
performance <- results |> 
  mutate(
    # compute lower and upper limits of a 95% confidence interval (ci_lower_XX, ci_upper_XX)
    # and check whether the population value lies within the 95% CI (covered_XX)
    # single group
    ci_lower_single_phi11_g1 = single_phi11_g1 - qnorm(0.975) * single_phi11_g1_se,
    ci_upper_single_phi11_g1 = single_phi11_g1 + qnorm(0.975) * single_phi11_g1_se,
    covered_single_phi11_g1 = (ci_lower_single_phi11_g1 <= phi11_g1_pop) & (ci_upper_single_phi11_g1 >= phi11_g1_pop),
    ci_lower_single_phi22_g1 = single_phi22_g1 - qnorm(0.975) * single_phi22_g1_se,
    ci_upper_single_phi22_g1 = single_phi22_g1 + qnorm(0.975) * single_phi22_g1_se,
    covered_single_phi22_g1 = (ci_lower_single_phi22_g1 <= phi22_g1_pop) & (ci_upper_single_phi22_g1 >= phi22_g1_pop),
    ci_lower_single_phi12_g1 = single_phi12_g1 - qnorm(0.975) * single_phi12_g1_se,
    ci_upper_single_phi12_g1 = single_phi12_g1 + qnorm(0.975) * single_phi12_g1_se,
    covered_single_phi12_g1 = (ci_lower_single_phi12_g1 <= phi12_g1_pop) & (ci_upper_single_phi12_g1 >= phi12_g1_pop),
    ci_lower_single_phi21_g1 = single_phi21_g1 - qnorm(0.975) * single_phi21_g1_se,
    ci_upper_single_phi21_g1 = single_phi21_g1 + qnorm(0.975) * single_phi21_g1_se,
    covered_single_phi21_g1 = (ci_lower_single_phi21_g1 <= phi21_g1_pop) & (ci_upper_single_phi21_g1 >= phi21_g1_pop),
    
    ci_lower_single_phi11_g2 = single_phi11_g2 - qnorm(0.975) * single_phi11_g2_se,
    ci_upper_single_phi11_g2 = single_phi11_g2 + qnorm(0.975) * single_phi11_g2_se,
    covered_single_phi11_g2 = (ci_lower_single_phi11_g2 <= phi11_g2_pop) & (ci_upper_single_phi11_g2 >= phi11_g2_pop),
    ci_lower_single_phi22_g2 = single_phi22_g2 - qnorm(0.975) * single_phi22_g2_se,
    ci_upper_single_phi22_g2 = single_phi22_g2 + qnorm(0.975) * single_phi22_g2_se,
    covered_single_phi22_g2 = (ci_lower_single_phi22_g2 <= phi22_g2_pop) & (ci_upper_single_phi22_g2 >= phi22_g2_pop),
    ci_lower_single_phi12_g2 = single_phi12_g2 - qnorm(0.975) * single_phi12_g2_se,
    ci_upper_single_phi12_g2 = single_phi12_g2 + qnorm(0.975) * single_phi12_g2_se,
    covered_single_phi12_g2 = (ci_lower_single_phi12_g2 <= phi12_g2_pop) & (ci_upper_single_phi12_g2 >= phi12_g2_pop),
    ci_lower_single_phi21_g2 = single_phi21_g2 - qnorm(0.975) * single_phi21_g2_se,
    ci_upper_single_phi21_g2 = single_phi21_g2 + qnorm(0.975) * single_phi21_g2_se,
    covered_single_phi21_g2 = (ci_lower_single_phi21_g2 <= phi21_g2_pop) & (ci_upper_single_phi21_g2 >= phi21_g2_pop),
    
    # multi group
    ci_lower_multi_phi11_g1 = multi_phi11_g1 - qnorm(0.975) * multi_phi11_g1_se,
    ci_upper_multi_phi11_g1 = multi_phi11_g1 + qnorm(0.975) * multi_phi11_g1_se,
    covered_multi_phi11_g1 = (ci_lower_multi_phi11_g1 <= phi11_g1_pop) & (ci_upper_multi_phi11_g1 >= phi11_g1_pop),
    ci_lower_multi_phi22_g1 = multi_phi22_g1 - qnorm(0.975) * multi_phi22_g1_se,
    ci_upper_multi_phi22_g1 = multi_phi22_g1 + qnorm(0.975) * multi_phi22_g1_se,
    covered_multi_phi22_g1 = (ci_lower_multi_phi22_g1 <= phi22_g1_pop) & (ci_upper_multi_phi22_g1 >= phi22_g1_pop),
    ci_lower_multi_phi12_g1 = multi_phi12_g1 - qnorm(0.975) * multi_phi12_g1_se,
    ci_upper_multi_phi12_g1 = multi_phi12_g1 + qnorm(0.975) * multi_phi12_g1_se,
    covered_multi_phi12_g1 = (ci_lower_multi_phi12_g1 <= phi12_g1_pop) & (ci_upper_multi_phi12_g1 >= phi12_g1_pop),
    ci_lower_multi_phi21_g1 = multi_phi21_g1 - qnorm(0.975) * multi_phi21_g1_se,
    ci_upper_multi_phi21_g1 = multi_phi21_g1 + qnorm(0.975) * multi_phi21_g1_se,
    covered_multi_phi21_g1 = (ci_lower_multi_phi21_g1 <= phi21_g1_pop) & (ci_upper_multi_phi21_g1 >= phi21_g1_pop),
    
    ci_lower_multi_phi11_g2 = multi_phi11_g2 - qnorm(0.975) * multi_phi11_g2_se,
    ci_upper_multi_phi11_g2 = multi_phi11_g2 + qnorm(0.975) * multi_phi11_g2_se,
    covered_multi_phi11_g2 = (ci_lower_multi_phi11_g2 <= phi11_g2_pop) & (ci_upper_multi_phi11_g2 >= phi11_g2_pop),
    ci_lower_multi_phi22_g2 = multi_phi22_g2 - qnorm(0.975) * multi_phi22_g2_se,
    ci_upper_multi_phi22_g2 = multi_phi22_g2 + qnorm(0.975) * multi_phi22_g2_se,
    covered_multi_phi22_g2 = (ci_lower_multi_phi22_g2 <= phi22_g2_pop) & (ci_upper_multi_phi22_g2 >= phi22_g2_pop),
    ci_lower_multi_phi12_g2 = multi_phi12_g2 - qnorm(0.975) * multi_phi12_g2_se,
    ci_upper_multi_phi12_g2 = multi_phi12_g2 + qnorm(0.975) * multi_phi12_g2_se,
    covered_multi_phi12_g2 = (ci_lower_multi_phi12_g2 <= phi12_g2_pop) & (ci_upper_multi_phi12_g2 >= phi12_g2_pop),
    ci_lower_multi_phi21_g2 = multi_phi21_g2 - qnorm(0.975) * multi_phi21_g2_se,
    ci_upper_multi_phi21_g2 = multi_phi21_g2 + qnorm(0.975) * multi_phi21_g2_se,
    covered_multi_phi21_g2 = (ci_lower_multi_phi21_g2 <= phi21_g2_pop) & (ci_upper_multi_phi21_g2 >= phi21_g2_pop),
  ) |> 
  
  # group by all simulation aspects (to obtain the unique conditions)
  group_by(invariance_level, pattern, ss_n, ss_t, ss_ratio) |> 
  
  # compute the outcomes per condition
  summarize(
    ## absolute bias (average estimate minus population value):
    AB_single_phi11_g1 = mean(single_phi11_g1) - unique(phi11_g1_pop),
    AB_single_phi22_g1 = mean(single_phi22_g1) - unique(phi22_g1_pop),
    AB_single_phi12_g1 = mean(single_phi12_g1) - unique(phi12_g1_pop),
    AB_single_phi21_g1 = mean(single_phi21_g1) - unique(phi21_g1_pop),
    
    AB_single_phi11_g2 = mean(single_phi11_g2) - unique(phi11_g2_pop),
    AB_single_phi22_g2 = mean(single_phi22_g2) - unique(phi22_g2_pop),
    AB_single_phi12_g2 = mean(single_phi12_g2) - unique(phi12_g2_pop),
    AB_single_phi21_g2 = mean(single_phi21_g2) - unique(phi21_g2_pop),
    
    AB_multi_phi11_g1 = mean(multi_phi11_g1) - unique(phi11_g1_pop),
    AB_multi_phi22_g1 = mean(multi_phi22_g1) - unique(phi22_g1_pop),
    AB_multi_phi12_g1 = mean(multi_phi12_g1) - unique(phi12_g1_pop),
    AB_multi_phi21_g1 = mean(multi_phi21_g1) - unique(phi21_g1_pop),
    
    AB_multi_phi11_g2 = mean(multi_phi11_g2) - unique(phi11_g2_pop),
    AB_multi_phi22_g2 = mean(multi_phi22_g2) - unique(phi22_g2_pop),
    AB_multi_phi12_g2 = mean(multi_phi12_g2) - unique(phi12_g2_pop),
    AB_multi_phi21_g2 = mean(multi_phi21_g2) - unique(phi21_g2_pop),
    
    ## RMSE
    RMSE_single_phi11_g1 = sqrt(mean((single_phi11_g1 - phi11_g1_pop)^2)),
    RMSE_single_phi22_g1 = sqrt(mean((single_phi22_g1 - phi22_g1_pop)^2)),
    RMSE_single_phi12_g1 = sqrt(mean((single_phi12_g1 - phi12_g1_pop)^2)),
    RMSE_single_phi21_g1 = sqrt(mean((single_phi21_g1 - phi21_g1_pop)^2)),
    
    RMSE_single_phi11_g2 = sqrt(mean((single_phi11_g2 - phi11_g2_pop)^2)),
    RMSE_single_phi22_g2 = sqrt(mean((single_phi22_g2 - phi22_g2_pop)^2)),
    RMSE_single_phi12_g2 = sqrt(mean((single_phi12_g2 - phi12_g2_pop)^2)),
    RMSE_single_phi21_g2 = sqrt(mean((single_phi21_g2 - phi21_g2_pop)^2)),
    
    RMSE_multi_phi11_g1 = sqrt(mean((multi_phi11_g1 - phi11_g1_pop)^2)),
    RMSE_multi_phi22_g1 = sqrt(mean((multi_phi22_g1 - phi22_g1_pop)^2)),
    RMSE_multi_phi12_g1 = sqrt(mean((multi_phi12_g1 - phi12_g1_pop)^2)),
    RMSE_multi_phi21_g1 = sqrt(mean((multi_phi21_g1 - phi21_g1_pop)^2)),
    
    RMSE_multi_phi11_g2 = sqrt(mean((multi_phi11_g2 - phi11_g2_pop)^2)),
    RMSE_multi_phi22_g2 = sqrt(mean((multi_phi22_g2 - phi22_g2_pop)^2)),
    RMSE_multi_phi12_g2 = sqrt(mean((multi_phi12_g2 - phi12_g2_pop)^2)),
    RMSE_multi_phi21_g2 = sqrt(mean((multi_phi21_g2 - phi21_g2_pop)^2)),
    
    ## Standard error recovery (average SE divided by empirical SE)
    SER_single_phi11_g1 = mean(single_phi11_g1_se) / sd(single_phi11_g1),
    SER_single_phi22_g1 = mean(single_phi22_g1_se) / sd(single_phi22_g1),
    SER_single_phi12_g1 = mean(single_phi12_g1_se) / sd(single_phi12_g1),
    SER_single_phi21_g1 = mean(single_phi21_g1_se) / sd(single_phi21_g1),
    
    SER_single_phi11_g2 = mean(single_phi11_g2_se) / sd(single_phi11_g2),
    SER_single_phi22_g2 = mean(single_phi22_g2_se) / sd(single_phi22_g2),
    SER_single_phi12_g2 = mean(single_phi12_g2_se) / sd(single_phi12_g2),
    SER_single_phi21_g2 = mean(single_phi21_g2_se) / sd(single_phi21_g2),
    
    SER_multi_phi11_g1 = mean(multi_phi11_g1_se) / sd(multi_phi11_g1),
    SER_multi_phi22_g1 = mean(multi_phi22_g1_se) / sd(multi_phi22_g1),
    SER_multi_phi12_g1 = mean(multi_phi12_g1_se) / sd(multi_phi12_g1),
    SER_multi_phi21_g1 = mean(multi_phi21_g1_se) / sd(multi_phi21_g1),
    
    SER_multi_phi11_g2 = mean(multi_phi11_g2_se) / sd(multi_phi11_g2),
    SER_multi_phi22_g2 = mean(multi_phi22_g2_se) / sd(multi_phi22_g2),
    SER_multi_phi12_g2 = mean(multi_phi12_g2_se) / sd(multi_phi12_g2),
    SER_multi_phi21_g2 = mean(multi_phi21_g2_se) / sd(multi_phi21_g2),
    
    ## coverage (proportion of replications where population value is within 95% CI)
    COV_single_phi11_g1 = mean(covered_single_phi11_g1),
    COV_single_phi22_g1 = mean(covered_single_phi22_g1),
    COV_single_phi12_g1 = mean(covered_single_phi12_g1),
    COV_single_phi21_g1 = mean(covered_single_phi21_g1),
    
    COV_single_phi11_g2 = mean(covered_single_phi11_g2),
    COV_single_phi22_g2 = mean(covered_single_phi22_g2),
    COV_single_phi12_g2 = mean(covered_single_phi12_g2),
    COV_single_phi21_g2 = mean(covered_single_phi21_g2),
    
    COV_multi_phi11_g1 = mean(covered_multi_phi11_g1),
    COV_multi_phi22_g1 = mean(covered_multi_phi22_g1),
    COV_multi_phi12_g1 = mean(covered_multi_phi12_g1),
    COV_multi_phi21_g1 = mean(covered_multi_phi21_g1),
    
    COV_multi_phi11_g2 = mean(covered_multi_phi11_g2),
    COV_multi_phi22_g2 = mean(covered_multi_phi22_g2),
    COV_multi_phi12_g2 = mean(covered_multi_phi12_g2),
    COV_multi_phi21_g2 = mean(covered_multi_phi21_g2),
    
    .groups = "drop") |> 
  
  # the following code transforms the columns
  # that outcomes and single/multi-group approach are now in different rows instead of columns
  # (introducing the column "type" = outcome, and "method" = modeling approach)
  pivot_longer(cols = 6:69,
               names_to = "key",
               values_to = "value") |> 
  separate_wider_delim(key, "_", names = c("type", "method", "parameter", "group")) |> 
  mutate(name = paste(parameter, group, sep = "_")) |> 
  dplyr::select(-parameter, -group) |> 
  pivot_wider(names_from = name,
              values_from = value)

performance_overall <- performance |> 
  group_by(method, type) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = NA,
         level = NA,
         .after = type)
print(performance_overall, width = Inf, n = Inf)

performance_invariance_level <- performance |> 
  group_by(method, type, invariance_level) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "measurement invariance",
         .after = type) |> 
  rename(level = invariance_level)
print(performance_invariance_level, width = Inf, n = Inf)

performance_pattern <- performance |> 
  group_by(method, type, pattern) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "pattern",
         .after = type) |> 
  rename(level = pattern)
print(performance_pattern, width = Inf, n = Inf)

performance_ss_n <- performance |> 
  group_by(method, type, ss_n) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "ss_n",
         .after = type) |> 
  rename(level = ss_n) |> 
  mutate(level = as.factor(level))
print(performance_ss_n, width = Inf, n = Inf)

performance_ss_t <- performance |> 
  group_by(method, type, ss_t) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "ss_t",
         .after = type) |> 
  rename(level = ss_t) |> 
  mutate(level = as.factor(level))
print(performance_ss_t, width = Inf, n = Inf)

performance_ss_ratio <- performance |> 
  group_by(method, type, ss_ratio) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "ss_ratio",
         .after = type) |> 
  rename(level = ss_ratio)
print(performance_ss_ratio, width = Inf, n = Inf)

performance_full <- bind_rows(performance_overall, performance_invariance_level, performance_pattern, performance_ss_n, performance_ss_t, performance_ss_ratio)
# write_csv(performance_full |> dplyr::filter(type == "AB", method == "multi"), file = "table_bias_multi.csv")
# write_csv(performance_full |> dplyr::filter(type == "RMSE", method == "multi"), file = "table_RMSE_multi.csv")
# write_csv(performance_full |> dplyr::filter(type == "COV", method == "multi"), file = "table_coverage_multi.csv")
# write_csv(performance_full |> dplyr::filter(type == "SER", method == "multi"), file = "table_SER_multi.csv")
# write_csv(performance_full |> dplyr::filter(type == "AB", method == "single"), file = "table_bias_single.csv")
# write_csv(performance_full |> dplyr::filter(type == "RMSE", method == "single"), file = "table_RMSE_single.csv")
# write_csv(performance_full |> dplyr::filter(type == "COV", method == "single"), file = "table_coverage_single.csv")
# write_csv(performance_full |> dplyr::filter(type == "SER", method == "single"), file = "table_SER_single.csv")


#### Simulation 2 ####
# read in the data
results <- read_csv("Simulation 2/sim2.csv",
                    col_types = cols(step1_personspec_warning_text = col_character(),
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
                                     step3_groupspec_error_text = col_character())) |> 
  mutate(heterogeneity = factor(heterogeneity, levels = c("small", "large")),
         nature = factor(nature, levels = c("continuous", "categorical"))
  ) |> 
  arrange(iteration)                                                            # sort by iteration

## Warnings and errors:
results |> 
  summarize(across(step1_personspec_warning:rerun_step1_groupspec, ~sum(.x, na.rm = TRUE))) |> 
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
    outlier_SE = if_any(c(
      personspec_phi11_g1_se, personspec_phi22_g1_se, personspec_phi12_g1_se, personspec_phi21_g1_se,
      personspec_phi11_g2_se, personspec_phi22_g2_se, personspec_phi12_g2_se, personspec_phi21_g2_se,
      groupspec_phi11_g1_se,  groupspec_phi22_g1_se,  groupspec_phi12_g1_se,  groupspec_phi21_g1_se,
      groupspec_phi11_g2_se,  groupspec_phi22_g2_se,  groupspec_phi12_g2_se,  groupspec_phi21_g2_se
    ), ~ .x > 1),
    # missing SEs (i.e., those that are NA)
    NA_SE = if_any(c(
      personspec_phi11_g1_se, personspec_phi22_g1_se, personspec_phi12_g1_se, personspec_phi21_g1_se,
      personspec_phi11_g2_se, personspec_phi22_g2_se, personspec_phi12_g2_se, personspec_phi21_g2_se,
      groupspec_phi11_g1_se,  groupspec_phi22_g1_se,  groupspec_phi12_g1_se,  groupspec_phi21_g1_se,
      groupspec_phi11_g2_se,  groupspec_phi22_g2_se,  groupspec_phi12_g2_se,  groupspec_phi21_g2_se
    ), ~ is.na(.x)))

sum(results$outlier_SE, na.rm = TRUE)
# 0 SEs > 1
sum(results$NA_SE)
# 7 missing SEs

results <- results |>
  mutate(remove = if_any(c(NA_SE, step1_personspec_warning, step3_groupspec_error)))
sum(results$remove, na.rm = TRUE)
# 20 data sets to remove

results <- results |>
  dplyr::filter(!remove)


# ## recovery of MM parameters:
# results |> 
#   summarize(mean_bias_lambda = mean(bias_lambda),
#             mean_bias_theta = mean(bias_theta),
#             mean_bias_tau = mean(bias_tau), 
#             mean_RMSE_lambda = mean(RMSE_lambda),
#             mean_RMSE_theta = mean(RMSE_theta),
#             mean_RMSE_tau = mean(RMSE_tau))

## recovery of SM parameters:
performance <- results |> 
  mutate(
    # compute lower and upper limits of a 95% confidence interval (ci_lower_XX, ci_upper_XX)
    # and check whether the population value lies within the 95% CI (covered_XX)
    # personspecific:
    ci_lower_personspec_phi11_g1 = personspec_phi11_g1 - qnorm(0.975) * personspec_phi11_g1_se,
    ci_upper_personspec_phi11_g1 = personspec_phi11_g1 + qnorm(0.975) * personspec_phi11_g1_se,
    covered_personspec_phi11_g1 = (ci_lower_personspec_phi11_g1 <= phi11_g1_pop) & (ci_upper_personspec_phi11_g1 >= phi11_g1_pop),
    ci_lower_personspec_phi22_g1 = personspec_phi22_g1 - qnorm(0.975) * personspec_phi22_g1_se,
    ci_upper_personspec_phi22_g1 = personspec_phi22_g1 + qnorm(0.975) * personspec_phi22_g1_se,
    covered_personspec_phi22_g1 = (ci_lower_personspec_phi22_g1 <= phi22_g1_pop) & (ci_upper_personspec_phi22_g1 >= phi22_g1_pop),
    ci_lower_personspec_phi12_g1 = personspec_phi12_g1 - qnorm(0.975) * personspec_phi12_g1_se,
    ci_upper_personspec_phi12_g1 = personspec_phi12_g1 + qnorm(0.975) * personspec_phi12_g1_se,
    covered_personspec_phi12_g1 = (ci_lower_personspec_phi12_g1 <= phi12_g1_pop) & (ci_upper_personspec_phi12_g1 >= phi12_g1_pop),
    ci_lower_personspec_phi21_g1 = personspec_phi21_g1 - qnorm(0.975) * personspec_phi21_g1_se,
    ci_upper_personspec_phi21_g1 = personspec_phi21_g1 + qnorm(0.975) * personspec_phi21_g1_se,
    covered_personspec_phi21_g1 = (ci_lower_personspec_phi21_g1 <= phi21_g1_pop) & (ci_upper_personspec_phi21_g1 >= phi21_g1_pop),
    
    ci_lower_personspec_phi11_g2 = personspec_phi11_g2 - qnorm(0.975) * personspec_phi11_g2_se,
    ci_upper_personspec_phi11_g2 = personspec_phi11_g2 + qnorm(0.975) * personspec_phi11_g2_se,
    covered_personspec_phi11_g2 = (ci_lower_personspec_phi11_g2 <= phi11_g2_pop) & (ci_upper_personspec_phi11_g2 >= phi11_g2_pop),
    ci_lower_personspec_phi22_g2 = personspec_phi22_g2 - qnorm(0.975) * personspec_phi22_g2_se,
    ci_upper_personspec_phi22_g2 = personspec_phi22_g2 + qnorm(0.975) * personspec_phi22_g2_se,
    covered_personspec_phi22_g2 = (ci_lower_personspec_phi22_g2 <= phi22_g2_pop) & (ci_upper_personspec_phi22_g2 >= phi22_g2_pop),
    ci_lower_personspec_phi12_g2 = personspec_phi12_g2 - qnorm(0.975) * personspec_phi12_g2_se,
    ci_upper_personspec_phi12_g2 = personspec_phi12_g2 + qnorm(0.975) * personspec_phi12_g2_se,
    covered_personspec_phi12_g2 = (ci_lower_personspec_phi12_g2 <= phi12_g2_pop) & (ci_upper_personspec_phi12_g2 >= phi12_g2_pop),
    ci_lower_personspec_phi21_g2 = personspec_phi21_g2 - qnorm(0.975) * personspec_phi21_g2_se,
    ci_upper_personspec_phi21_g2 = personspec_phi21_g2 + qnorm(0.975) * personspec_phi21_g2_se,
    covered_personspec_phi21_g2 = (ci_lower_personspec_phi21_g2 <= phi21_g2_pop) & (ci_upper_personspec_phi21_g2 >= phi21_g2_pop),
    
    # group-specific
    ci_lower_groupspec_phi11_g1 = groupspec_phi11_g1 - qnorm(0.975) * groupspec_phi11_g1_se,
    ci_upper_groupspec_phi11_g1 = groupspec_phi11_g1 + qnorm(0.975) * groupspec_phi11_g1_se,
    covered_groupspec_phi11_g1 = (ci_lower_groupspec_phi11_g1 <= phi11_g1_pop) & (ci_upper_groupspec_phi11_g1 >= phi11_g1_pop),
    ci_lower_groupspec_phi22_g1 = groupspec_phi22_g1 - qnorm(0.975) * groupspec_phi22_g1_se,
    ci_upper_groupspec_phi22_g1 = groupspec_phi22_g1 + qnorm(0.975) * groupspec_phi22_g1_se,
    covered_groupspec_phi22_g1 = (ci_lower_groupspec_phi22_g1 <= phi22_g1_pop) & (ci_upper_groupspec_phi22_g1 >= phi22_g1_pop),
    ci_lower_groupspec_phi12_g1 = groupspec_phi12_g1 - qnorm(0.975) * groupspec_phi12_g1_se,
    ci_upper_groupspec_phi12_g1 = groupspec_phi12_g1 + qnorm(0.975) * groupspec_phi12_g1_se,
    covered_groupspec_phi12_g1 = (ci_lower_groupspec_phi12_g1 <= phi12_g1_pop) & (ci_upper_groupspec_phi12_g1 >= phi12_g1_pop),
    ci_lower_groupspec_phi21_g1 = groupspec_phi21_g1 - qnorm(0.975) * groupspec_phi21_g1_se,
    ci_upper_groupspec_phi21_g1 = groupspec_phi21_g1 + qnorm(0.975) * groupspec_phi21_g1_se,
    covered_groupspec_phi21_g1 = (ci_lower_groupspec_phi21_g1 <= phi21_g1_pop) & (ci_upper_groupspec_phi21_g1 >= phi21_g1_pop),
    
    ci_lower_groupspec_phi11_g2 = groupspec_phi11_g2 - qnorm(0.975) * groupspec_phi11_g2_se,
    ci_upper_groupspec_phi11_g2 = groupspec_phi11_g2 + qnorm(0.975) * groupspec_phi11_g2_se,
    covered_groupspec_phi11_g2 = (ci_lower_groupspec_phi11_g2 <= phi11_g2_pop) & (ci_upper_groupspec_phi11_g2 >= phi11_g2_pop),
    ci_lower_groupspec_phi22_g2 = groupspec_phi22_g2 - qnorm(0.975) * groupspec_phi22_g2_se,
    ci_upper_groupspec_phi22_g2 = groupspec_phi22_g2 + qnorm(0.975) * groupspec_phi22_g2_se,
    covered_groupspec_phi22_g2 = (ci_lower_groupspec_phi22_g2 <= phi22_g2_pop) & (ci_upper_groupspec_phi22_g2 >= phi22_g2_pop),
    ci_lower_groupspec_phi12_g2 = groupspec_phi12_g2 - qnorm(0.975) * groupspec_phi12_g2_se,
    ci_upper_groupspec_phi12_g2 = groupspec_phi12_g2 + qnorm(0.975) * groupspec_phi12_g2_se,
    covered_groupspec_phi12_g2 = (ci_lower_groupspec_phi12_g2 <= phi12_g2_pop) & (ci_upper_groupspec_phi12_g2 >= phi12_g2_pop),
    ci_lower_groupspec_phi21_g2 = groupspec_phi21_g2 - qnorm(0.975) * groupspec_phi21_g2_se,
    ci_upper_groupspec_phi21_g2 = groupspec_phi21_g2 + qnorm(0.975) * groupspec_phi21_g2_se,
    covered_groupspec_phi21_g2 = (ci_lower_groupspec_phi21_g2 <= phi21_g2_pop) & (ci_upper_groupspec_phi21_g2 >= phi21_g2_pop),
  ) |> 
  
  # group by all simulation aspects (to obtain the unique conditions)
  group_by(ss_t, heterogeneity, nature) |> 
  
  # compute the outcomes per condition
  summarize(
    ## absolute bias (average estimate minus population value):
    AB_personspec_phi11_g1 = mean(personspec_phi11_g1) - unique(phi11_g1_pop),
    AB_personspec_phi22_g1 = mean(personspec_phi22_g1) - unique(phi22_g1_pop),
    AB_personspec_phi12_g1 = mean(personspec_phi12_g1) - unique(phi12_g1_pop),
    AB_personspec_phi21_g1 = mean(personspec_phi21_g1) - unique(phi21_g1_pop),
    
    AB_personspec_phi11_g2 = mean(personspec_phi11_g2) - unique(phi11_g2_pop),
    AB_personspec_phi22_g2 = mean(personspec_phi22_g2) - unique(phi22_g2_pop),
    AB_personspec_phi12_g2 = mean(personspec_phi12_g2) - unique(phi12_g2_pop),
    AB_personspec_phi21_g2 = mean(personspec_phi21_g2) - unique(phi21_g2_pop),
    
    AB_groupspec_phi11_g1 = mean(groupspec_phi11_g1) - unique(phi11_g1_pop),
    AB_groupspec_phi22_g1 = mean(groupspec_phi22_g1) - unique(phi22_g1_pop),
    AB_groupspec_phi12_g1 = mean(groupspec_phi12_g1) - unique(phi12_g1_pop),
    AB_groupspec_phi21_g1 = mean(groupspec_phi21_g1) - unique(phi21_g1_pop),
    
    AB_groupspec_phi11_g2 = mean(groupspec_phi11_g2) - unique(phi11_g2_pop),
    AB_groupspec_phi22_g2 = mean(groupspec_phi22_g2) - unique(phi22_g2_pop),
    AB_groupspec_phi12_g2 = mean(groupspec_phi12_g2) - unique(phi12_g2_pop),
    AB_groupspec_phi21_g2 = mean(groupspec_phi21_g2) - unique(phi21_g2_pop),
    
    ## RMSE
    RMSE_personspec_phi11_g1 = sqrt(mean((personspec_phi11_g1 - phi11_g1_pop)^2)),
    RMSE_personspec_phi22_g1 = sqrt(mean((personspec_phi22_g1 - phi22_g1_pop)^2)),
    RMSE_personspec_phi12_g1 = sqrt(mean((personspec_phi12_g1 - phi12_g1_pop)^2)),
    RMSE_personspec_phi21_g1 = sqrt(mean((personspec_phi21_g1 - phi21_g1_pop)^2)),
    
    RMSE_personspec_phi11_g2 = sqrt(mean((personspec_phi11_g2 - phi11_g2_pop)^2)),
    RMSE_personspec_phi22_g2 = sqrt(mean((personspec_phi22_g2 - phi22_g2_pop)^2)),
    RMSE_personspec_phi12_g2 = sqrt(mean((personspec_phi12_g2 - phi12_g2_pop)^2)),
    RMSE_personspec_phi21_g2 = sqrt(mean((personspec_phi21_g2 - phi21_g2_pop)^2)),
    
    RMSE_groupspec_phi11_g1 = sqrt(mean((groupspec_phi11_g1 - phi11_g1_pop)^2)),
    RMSE_groupspec_phi22_g1 = sqrt(mean((groupspec_phi22_g1 - phi22_g1_pop)^2)),
    RMSE_groupspec_phi12_g1 = sqrt(mean((groupspec_phi12_g1 - phi12_g1_pop)^2)),
    RMSE_groupspec_phi21_g1 = sqrt(mean((groupspec_phi21_g1 - phi21_g1_pop)^2)),
    
    RMSE_groupspec_phi11_g2 = sqrt(mean((groupspec_phi11_g2 - phi11_g2_pop)^2)),
    RMSE_groupspec_phi22_g2 = sqrt(mean((groupspec_phi22_g2 - phi22_g2_pop)^2)),
    RMSE_groupspec_phi12_g2 = sqrt(mean((groupspec_phi12_g2 - phi12_g2_pop)^2)),
    RMSE_groupspec_phi21_g2 = sqrt(mean((groupspec_phi21_g2 - phi21_g2_pop)^2)),
    
    ## Standard error recovery (average SE divided by empirical SE)
    SER_personspec_phi11_g1 = mean(personspec_phi11_g1_se) / sd(personspec_phi11_g1),
    SER_personspec_phi22_g1 = mean(personspec_phi22_g1_se) / sd(personspec_phi22_g1),
    SER_personspec_phi12_g1 = mean(personspec_phi12_g1_se) / sd(personspec_phi12_g1),
    SER_personspec_phi21_g1 = mean(personspec_phi21_g1_se) / sd(personspec_phi21_g1),
    
    SER_personspec_phi11_g2 = mean(personspec_phi11_g2_se) / sd(personspec_phi11_g2),
    SER_personspec_phi22_g2 = mean(personspec_phi22_g2_se) / sd(personspec_phi22_g2),
    SER_personspec_phi12_g2 = mean(personspec_phi12_g2_se) / sd(personspec_phi12_g2),
    SER_personspec_phi21_g2 = mean(personspec_phi21_g2_se) / sd(personspec_phi21_g2),
    
    SER_groupspec_phi11_g1 = mean(groupspec_phi11_g1_se) / sd(groupspec_phi11_g1),
    SER_groupspec_phi22_g1 = mean(groupspec_phi22_g1_se) / sd(groupspec_phi22_g1),
    SER_groupspec_phi12_g1 = mean(groupspec_phi12_g1_se) / sd(groupspec_phi12_g1),
    SER_groupspec_phi21_g1 = mean(groupspec_phi21_g1_se) / sd(groupspec_phi21_g1),
    
    SER_groupspec_phi11_g2 = mean(groupspec_phi11_g2_se) / sd(groupspec_phi11_g2),
    SER_groupspec_phi22_g2 = mean(groupspec_phi22_g2_se) / sd(groupspec_phi22_g2),
    SER_groupspec_phi12_g2 = mean(groupspec_phi12_g2_se) / sd(groupspec_phi12_g2),
    SER_groupspec_phi21_g2 = mean(groupspec_phi21_g2_se) / sd(groupspec_phi21_g2),
    
    ## coverage (proportion of replications where population value is within 95% CI)
    COV_personspec_phi11_g1 = mean(covered_personspec_phi11_g1),
    COV_personspec_phi22_g1 = mean(covered_personspec_phi22_g1),
    COV_personspec_phi12_g1 = mean(covered_personspec_phi12_g1),
    COV_personspec_phi21_g1 = mean(covered_personspec_phi21_g1),
    
    COV_personspec_phi11_g2 = mean(covered_personspec_phi11_g2),
    COV_personspec_phi22_g2 = mean(covered_personspec_phi22_g2),
    COV_personspec_phi12_g2 = mean(covered_personspec_phi12_g2),
    COV_personspec_phi21_g2 = mean(covered_personspec_phi21_g2),
    
    COV_groupspec_phi11_g1 = mean(covered_groupspec_phi11_g1),
    COV_groupspec_phi22_g1 = mean(covered_groupspec_phi22_g1),
    COV_groupspec_phi12_g1 = mean(covered_groupspec_phi12_g1),
    COV_groupspec_phi21_g1 = mean(covered_groupspec_phi21_g1),
    
    COV_groupspec_phi11_g2 = mean(covered_groupspec_phi11_g2),
    COV_groupspec_phi22_g2 = mean(covered_groupspec_phi22_g2),
    COV_groupspec_phi12_g2 = mean(covered_groupspec_phi12_g2),
    COV_groupspec_phi21_g2 = mean(covered_groupspec_phi21_g2),
    
    .groups = "drop") |> 
  
  # the following code transforms the columns
  # that outcomes and single/multi-group approach are now in different rows instead of columns
  # (introducing the column "type" = outcome, and "method" = modeling approach)
  pivot_longer(cols = 4:67,
               names_to = "key",
               values_to = "value") |> 
  separate_wider_delim(key, "_", names = c("type", "method", "parameter", "group")) |> 
  mutate(name = paste(parameter, group, sep = "_")) |> 
  dplyr::select(-parameter, -group) |> 
  pivot_wider(names_from = name,
              values_from = value)

performance_overall <- performance |> 
  group_by(method, type) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = NA,
         level = NA,
         .after = type)
print(performance_overall, width = Inf, n = Inf)

performance_ss_t <- performance |> 
  group_by(method, type, ss_t) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "ss_t",
         .after = type) |> 
  rename(level = ss_t) |> 
  mutate(level = as.factor(level))
print(performance_ss_t, width = Inf, n = Inf)

performance_heterogeneity <- performance |> 
  group_by(method, type, heterogeneity) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "heterogeneity",
         .after = type) |> 
  rename(level = heterogeneity)
print(performance_heterogeneity, width = Inf, n = Inf)

performance_nature <- performance |> 
  group_by(method, type, nature) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |> 
  mutate(manipulated_aspect = "nature",
         .after = type) |> 
  rename(level = nature)
print(performance_nature, width = Inf, n = Inf)

performance_full <- bind_rows(performance_overall, performance_ss_t, performance_heterogeneity, performance_nature)
# write_csv(performance_full |> dplyr::filter(type == "AB", method == "personspec"), file = "table_bias_personspec.csv")
# write_csv(performance_full |> dplyr::filter(type == "RMSE", method == "personspec"), file = "table_RMSE_personspec.csv")
# write_csv(performance_full |> dplyr::filter(type == "COV", method == "personspec"), file = "table_coverage_personspec.csv")
# write_csv(performance_full |> dplyr::filter(type == "SER", method == "personspec"), file = "table_SER_personspec.csv")
# write_csv(performance_full |> dplyr::filter(type == "AB", method == "groupspec"), file = "table_bias_groupspec.csv")
# write_csv(performance_full |> dplyr::filter(type == "RMSE", method == "groupspec"), file = "table_RMSE_groupspec.csv")
# write_csv(performance_full |> dplyr::filter(type == "COV", method == "groupspec"), file = "table_coverage_groupspec.csv")
# write_csv(performance_full |> dplyr::filter(type == "SER", method == "groupspec"), file = "table_SER_groupspec.csv")