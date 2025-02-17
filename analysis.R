#### This script performs the analysis on the simulation output ####
library(tidyverse)
results <- read_csv("sim1.csv",
                    col_types = cols(step1_warning_text = col_character(),
                                     step2_warning_text = col_character(),
                                     step3_warning_text = col_character(),
                                     step1_error_text = col_character(),
                                     step2_error_text = col_character(),
                                     step3_error_text = col_character())) |> 
  mutate(lambda_noninvariance = factor(lambda_noninvariance, levels = c("no", "uniform", "mixed")),
         theta_noninvariance = factor(theta_noninvariance, levels = c("no", "uniform", "mixed")),
         nu_noninvariance = factor(nu_noninvariance, levels = c("no", "uniform", "mixed")),
         step1modeling = factor(step1modeling, levels = c("single-group", "multi-group"))
  ) %>% 
  arrange(iteration)                                                            # sort by iteration

results |> 
  group_by(step1modeling) |> 
  summarize(across(step1_warning:step3_error, ~sum(.x, na.rm = TRUE)))

unique(results$step3_warning_text)

results <- results |> 
  dplyr::filter(!step1_warning, !step2_warning, !step3_warning,
                !step1_error, !step2_error, !step3_error)

performance <- results |> 
  group_by(n, obs, lambda_noninvariance, theta_noninvariance, nu_noninvariance, step1modeling) |> 
  summarize(AB_phi11_g1 = (mean(phi11_g1, na.rm = TRUE) - unique(phi11_g1_pop))/unique(phi11_g1_pop),
            AB_phi22_g1 = (mean(phi22_g1, na.rm = TRUE) - unique(phi22_g1_pop))/unique(phi22_g1_pop),
            AB_phi12_g1 = (mean(phi12_g1, na.rm = TRUE) - unique(phi12_g1_pop))/unique(phi12_g1_pop),
            AB_phi21_g1 = (mean(phi21_g1, na.rm = TRUE) - unique(phi21_g1_pop))/unique(phi21_g1_pop),
            
            AB_phi11_g2 = (mean(phi11_g2, na.rm = TRUE) - unique(phi11_g2_pop))/unique(phi11_g2_pop),
            AB_phi22_g2 = (mean(phi22_g2, na.rm = TRUE) - unique(phi22_g2_pop))/unique(phi22_g2_pop),
            AB_phi12_g2 = (mean(phi12_g2, na.rm = TRUE) - unique(phi12_g2_pop))/unique(phi12_g2_pop),
            AB_phi21_g2 = (mean(phi21_g2, na.rm = TRUE) - unique(phi21_g2_pop))/unique(phi21_g2_pop),
            .groups = "drop")

performance_overall <- performance |> 
  group_by(step1modeling) |> 
  summarise(across(AB_phi11_g1:AB_phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_overall)

performance_n <- performance |> 
  group_by(n, step1modeling) |> 
  summarise(across(AB_phi11_g1:AB_phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_n)

performance_obs <- performance |> 
  group_by(obs, step1modeling) |> 
  summarise(across(AB_phi11_g1:AB_phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_obs)

performance_lambda_noninvariance <- performance |> 
  group_by(lambda_noninvariance, step1modeling) |> 
  summarise(across(AB_phi11_g1:AB_phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_lambda_noninvariance)

performance_theta_noninvariance <- performance |> 
  group_by(theta_noninvariance, step1modeling) |> 
  summarise(across(AB_phi11_g1:AB_phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_theta_noninvariance)

performance_nu_noninvariance <- performance |> 
  group_by(nu_noninvariance,step1modeling) |> 
  summarise(across(AB_phi11_g1:AB_phi21_g2, ~ mean(.x, na.rm = TRUE)))
print(performance_nu_noninvariance)



