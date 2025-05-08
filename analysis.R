#### This script performs the analysis on the simulation output ####
library(tidyverse)

#### Simulation 1 ####
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

## Warnings and errors:
results |> 
  summarize(across(step1_single_warning:step3_multi_error, ~sum(.x, na.rm = TRUE))) |> 
  print(width = Inf)
# no warnings or errors

results <- results |> 
  filter(!step1_single_warning)

results |> group_by(ss_n, ss_t) |> 
  summarize(mean_bias_lambda = mean(bias_lambda),
            mean_bias_theta = mean(bias_theta),
            mean_bias_tau = mean(bias_tau),
            mean_RMSE_lambda = mean(RMSE_lambda),
            mean_RMSE_theta = mean(RMSE_theta),
            mean_RMSE_tau = mean(RMSE_tau))

performance <- results |> 
  group_by(invariance_level, pattern, ss_n, ss_t, ss_ratio) |> 
  summarize(AB_single_phi11_g1 = (mean(single_phi11_g1, na.rm = TRUE) - unique(phi11_g1_pop)),
            AB_single_phi22_g1 = (mean(single_phi22_g1, na.rm = TRUE) - unique(phi22_g1_pop)),
            AB_single_phi12_g1 = (mean(single_phi12_g1, na.rm = TRUE) - unique(phi12_g1_pop)),
            AB_single_phi21_g1 = (mean(single_phi21_g1, na.rm = TRUE) - unique(phi21_g1_pop)),
            
            AB_single_phi11_g2 = (mean(single_phi11_g2, na.rm = TRUE) - unique(phi11_g2_pop)),
            AB_single_phi22_g2 = (mean(single_phi22_g2, na.rm = TRUE) - unique(phi22_g2_pop)),
            AB_single_phi12_g2 = (mean(single_phi12_g2, na.rm = TRUE) - unique(phi12_g2_pop)),
            AB_single_phi21_g2 = (mean(single_phi21_g2, na.rm = TRUE) - unique(phi21_g2_pop)),
            
            AB_multi_phi11_g1 = (mean(multi_phi11_g1, na.rm = TRUE) - unique(phi11_g1_pop)),
            AB_multi_phi22_g1 = (mean(multi_phi22_g1, na.rm = TRUE) - unique(phi22_g1_pop)),
            AB_multi_phi12_g1 = (mean(multi_phi12_g1, na.rm = TRUE) - unique(phi12_g1_pop)),
            AB_multi_phi21_g1 = (mean(multi_phi21_g1, na.rm = TRUE) - unique(phi21_g1_pop)),
            
            AB_multi_phi11_g2 = (mean(multi_phi11_g2, na.rm = TRUE) - unique(phi11_g2_pop)),
            AB_multi_phi22_g2 = (mean(multi_phi22_g2, na.rm = TRUE) - unique(phi22_g2_pop)),
            AB_multi_phi12_g2 = (mean(multi_phi12_g2, na.rm = TRUE) - unique(phi12_g2_pop)),
            AB_multi_phi21_g2 = (mean(multi_phi21_g2, na.rm = TRUE) - unique(phi21_g2_pop)),
            
            RB_single_phi11_g1 = (mean(single_phi11_g1, na.rm = TRUE) - unique(phi11_g1_pop))/unique(phi11_g1_pop),
            RB_single_phi22_g1 = (mean(single_phi22_g1, na.rm = TRUE) - unique(phi22_g1_pop))/unique(phi22_g1_pop),
            RB_single_phi12_g1 = (mean(single_phi12_g1, na.rm = TRUE) - unique(phi12_g1_pop))/unique(phi12_g1_pop),
            RB_single_phi21_g1 = (mean(single_phi21_g1, na.rm = TRUE) - unique(phi21_g1_pop))/unique(phi21_g1_pop),
            
            RB_single_phi11_g2 = (mean(single_phi11_g2, na.rm = TRUE) - unique(phi11_g2_pop))/unique(phi11_g2_pop),
            RB_single_phi22_g2 = (mean(single_phi22_g2, na.rm = TRUE) - unique(phi22_g2_pop))/unique(phi22_g2_pop),
            RB_single_phi12_g2 = (mean(single_phi12_g2, na.rm = TRUE) - unique(phi12_g2_pop))/unique(phi12_g2_pop),
            RB_single_phi21_g2 = (mean(single_phi21_g2, na.rm = TRUE) - unique(phi21_g2_pop))/unique(phi21_g2_pop),
            
            RB_multi_phi11_g1 = (mean(multi_phi11_g1, na.rm = TRUE) - unique(phi11_g1_pop))/unique(phi11_g1_pop),
            RB_multi_phi22_g1 = (mean(multi_phi22_g1, na.rm = TRUE) - unique(phi22_g1_pop))/unique(phi22_g1_pop),
            RB_multi_phi12_g1 = (mean(multi_phi12_g1, na.rm = TRUE) - unique(phi12_g1_pop))/unique(phi12_g1_pop),
            RB_multi_phi21_g1 = (mean(multi_phi21_g1, na.rm = TRUE) - unique(phi21_g1_pop))/unique(phi21_g1_pop),
            
            RB_multi_phi11_g2 = (mean(multi_phi11_g2, na.rm = TRUE) - unique(phi11_g2_pop))/unique(phi11_g2_pop),
            RB_multi_phi22_g2 = (mean(multi_phi22_g2, na.rm = TRUE) - unique(phi22_g2_pop))/unique(phi22_g2_pop),
            RB_multi_phi12_g2 = (mean(multi_phi12_g2, na.rm = TRUE) - unique(phi12_g2_pop))/unique(phi12_g2_pop),
            RB_multi_phi21_g2 = (mean(multi_phi21_g2, na.rm = TRUE) - unique(phi21_g2_pop))/unique(phi21_g2_pop),
            .groups = "drop") |> 
  pivot_longer(cols = 6:37,
               names_to = "key",
               values_to = "value") |> 
  separate_wider_delim(key, "_", names = c("type", "method", "parameter", "group")) |> 
  mutate(name = paste(parameter, group, sep = "_")) |> 
  dplyr::select(-parameter, -group) |> 
  pivot_wider(names_from = name,
              values_from = value)

performance_overall <- performance |> 
  group_by(type, method) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_overall, width = Inf, n = Inf)

performance_invariance_level <- performance |> 
  group_by(type, method, invariance_level) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_invariance_level, width = Inf, n = Inf)

performance_pattern <- performance |> 
  group_by(type, method, pattern) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_pattern, width = Inf, n = Inf)

performance_ss_n <- performance |> 
  group_by(type, method, invariance_level, ss_n) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_ss_n, width = Inf, n = Inf)

performance_ss_t <- performance |> 
  group_by(type, method, invariance_level, ss_t) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_ss_t, width = Inf, n = Inf)

performance_ss_ratio <- performance |> 
  group_by(type, method, ss_ratio) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_ss_ratio, width = Inf, n = Inf)
