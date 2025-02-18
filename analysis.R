#### This script performs the analysis on the simulation output ####
library(tidyverse)
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
  mutate(lambda_noninvariance = factor(lambda_noninvariance, levels = c("no", "uniform", "mixed")),
         theta_noninvariance = factor(theta_noninvariance, levels = c("no", "uniform", "mixed")),
         tau_noninvariance = factor(tau_noninvariance, levels = c("no", "uniform", "mixed"))
  ) |> 
  arrange(iteration)                                                            # sort by iteration

results |> 
  summarize(across(step1_single_warning:step3_multi_error, ~sum(.x, na.rm = TRUE))) |> 
  print(width = Inf)

unique(c(results$step3_single_error_text, results$step3_multi_warning_text))

results$pos[which(results$step3_single_warning)]

results <- results |> 
  dplyr::filter(!step3_single_warning, !step3_multi_warning)

performance <- results |> 
  group_by(n, obs, lambda_noninvariance, theta_noninvariance, tau_noninvariance) |> 
  summarize(AB_single_phi11_g1 = (mean(single_phi11_g1, na.rm = TRUE) - unique(phi11_g1_pop))/unique(phi11_g1_pop),
            AB_single_phi22_g1 = (mean(single_phi22_g1, na.rm = TRUE) - unique(phi22_g1_pop))/unique(phi22_g1_pop),
            AB_single_phi12_g1 = (mean(single_phi12_g1, na.rm = TRUE) - unique(phi12_g1_pop))/unique(phi12_g1_pop),
            AB_single_phi21_g1 = (mean(single_phi21_g1, na.rm = TRUE) - unique(phi21_g1_pop))/unique(phi21_g1_pop),
            
            AB_single_phi11_g2 = (mean(single_phi11_g2, na.rm = TRUE) - unique(phi11_g2_pop))/unique(phi11_g2_pop),
            AB_single_phi22_g2 = (mean(single_phi22_g2, na.rm = TRUE) - unique(phi22_g2_pop))/unique(phi22_g2_pop),
            AB_single_phi12_g2 = (mean(single_phi12_g2, na.rm = TRUE) - unique(phi12_g2_pop))/unique(phi12_g2_pop),
            AB_single_phi21_g2 = (mean(single_phi21_g2, na.rm = TRUE) - unique(phi21_g2_pop))/unique(phi21_g2_pop),
            
            AB_multi_phi11_g1 = (mean(multi_phi11_g1, na.rm = TRUE) - unique(phi11_g1_pop))/unique(phi11_g1_pop),
            AB_multi_phi22_g1 = (mean(multi_phi22_g1, na.rm = TRUE) - unique(phi22_g1_pop))/unique(phi22_g1_pop),
            AB_multi_phi12_g1 = (mean(multi_phi12_g1, na.rm = TRUE) - unique(phi12_g1_pop))/unique(phi12_g1_pop),
            AB_multi_phi21_g1 = (mean(multi_phi21_g1, na.rm = TRUE) - unique(phi21_g1_pop))/unique(phi21_g1_pop),
            
            AB_multi_phi11_g2 = (mean(multi_phi11_g2, na.rm = TRUE) - unique(phi11_g2_pop))/unique(phi11_g2_pop),
            AB_multi_phi22_g2 = (mean(multi_phi22_g2, na.rm = TRUE) - unique(phi22_g2_pop))/unique(phi22_g2_pop),
            AB_multi_phi12_g2 = (mean(multi_phi12_g2, na.rm = TRUE) - unique(phi12_g2_pop))/unique(phi12_g2_pop),
            AB_multi_phi21_g2 = (mean(multi_phi21_g2, na.rm = TRUE) - unique(phi21_g2_pop))/unique(phi21_g2_pop),
            .groups = "drop") |> 
  pivot_longer(cols = 6:21,
               names_to = "key",
               values_to = "value") |> 
  separate_wider_delim(key, "_", names = c("type", "method", "parameter", "group")) |> 
  mutate(name = paste(parameter, group, sep = "_")) |> 
  dplyr::select(-parameter, -group) |> 
  pivot_wider(names_from = name,
              values_from = value)

performance_overall <- performance |> 
  group_by(method) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_overall, width = Inf)

performance_n <- performance |> 
  group_by(method, n) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_n, width = Inf)

performance_obs <- performance |> 
  group_by(method, obs) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_obs, width = Inf)

performance_lambda_noninvariance <- performance |> 
  group_by(method, lambda_noninvariance) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_lambda_noninvariance, width = Inf)

performance_theta_noninvariance <- performance |> 
  group_by(method, theta_noninvariance) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_theta_noninvariance, width = Inf)

performance_tau_noninvariance <- performance |> 
  group_by(method, tau_noninvariance) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_tau_noninvariance, width = Inf)


performance_all_invariances <- performance |> 
  group_by(lambda_noninvariance, theta_noninvariance, tau_noninvariance, method) |> 
  summarise(across(phi11_g1:phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_all_invariances, width = Inf)

write_csv(performance_all_invariances, file = "test.csv")
