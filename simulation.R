#### This script performs the simulation and saves the output ####

#### load packages and functions ####
library(conflicted)
library(dplyr)
library(flock)
library(lavaan)
library(MASS)
library(OpenMx)
library(purrr)
library(tidyr)
library(stringr)

# packages for parallel processing and data reading:
library(parabar)
library(parallel)
library(readr)

# load functions
source("step1.R")
source("step2.R")
source("step3.R")
source("auxiliary_functions.R")


#### Define condition grids ####
## Simulation 1:
cond_sim1 <- expand.grid(replication = 1:200,
                    invariance_level = c("full_strict", "full_scalar", "partial_scalar", "partial_metric"),
                    pattern = c("unidirectional", "mixed"),
                    ss_n = c(24, 72),
                    ss_t = c(14, 56),
                    ss_ratio = c("balanced", "unbalanced"))


# add seeds:
set.seed(123)
cond_sim1$seed <- sample(1:(nrow(cond_sim1)*5), size = nrow(cond_sim1), replace = FALSE)
# add iteration number:
cond_sim1$iteration <- 1:nrow(cond_sim1)                                                  # (unique) number of each iteration


## Simulation 2:
cond_sim2 <- expand.grid(replication = 1:200,
                         ss_t = c(56, 84, 112),
                         heterogeneity = c("small", "large"),
                         nature = c("categorical", "continuous"))


# add seeds:
set.seed(456)
cond_sim2$seed <- sample(1:(nrow(cond_sim2)*5), size = nrow(cond_sim2), replace = FALSE)
# add iteration number:
cond_sim2$iteration <- 1:nrow(cond_sim2)                                                  # (unique) number of each iteration

#### Simulation Study 1 ####
source("Simulation 1/do_sim1.R")

## set up parallel computing
# open cluster
numCores <- parallel::detectCores() - 1
backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")


# load libraries in cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(flock)
  library(lavaan)
  library(MASS)
  library(OpenMx)
  library(purrr)
  library(tidyr)
  library(stringr)
  
  source("step1.R")
  source("step2.R")
  source("step3.R")
  source("Simulation 1/do_sim1.R")
  source("auxiliary_functions.R")
})

## load objects in cluster
export(backend, "cond_sim1")

## perform simulation
start  <- Sys.time()
output <- par_lapply(backend, 1:nrow(cond_sim1), do_sim1, cond = cond_sim1,
                     outputfile = "Simulation 1/sim1.csv", verbose = FALSE)
end <- Sys.time()

# close cluster:
stop_backend(backend)

#### Simulation Study 1 - Reestimation ####
source("Simulation 1/do_sim1_reestimation.R")

# load data set
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

sum(results$reestimate)
# there are in total 1183 data sets where there was either a warning (negative variances), error (failed step 3 estimation), or SE that was NA or greater than 1
table(results$reestimate, results$invariance_level)
table(results$reestimate, results$pattern)
table(results$reestimate, results$ss_n)
table(results$reestimate, results$ss_t)
table(results$reestimate, results$ss_ratio)

# get iteration numbers that need to be reestimated
new_iterations <- results$iteration[results$reestimate]

cond_sim1_reestimation <- cond_sim1[cond_sim1$iteration %in% new_iterations,]


## set up parallel computing
# open cluster
numCores <- parallel::detectCores() - 1
backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")


# load libraries in cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(flock)
  library(lavaan)
  library(MASS)
  library(OpenMx)
  library(purrr)
  library(tidyr)
  library(stringr)
  
  source("step1.R")
  source("step2.R")
  source("step3.R")
  source("Simulation 1/do_sim1_reestimation.R")
  source("auxiliary_functions.R")
})

## load objects in cluster
export(backend, "cond_sim1_reestimation")

## perform simulation
start  <- Sys.time()
output <- par_lapply(backend, 1:nrow(cond_sim1_reestimation), do_sim1_reestimation, cond = cond_sim1_reestimation,
                     outputfile = "Simulation 1/sim1_reestimation.csv", verbose = FALSE)
end <- Sys.time()

# close cluster:
stop_backend(backend)

#### Simulation Study 2 ####
source("Simulation 2/do_sim2.R")

## set up parallel computing
# open cluster
numCores <- parallel::detectCores() - 1
backend <- start_backend(numCores, cluster_type = "psock", backend_type = "async")


# load libraries in cluster
parabar::evaluate(backend, {
  library(conflicted)
  library(dplyr)
  library(flock)
  library(lavaan)
  library(MASS)
  library(OpenMx)
  library(purrr)
  library(tidyr)
  library(stringr)
  
  source("step1.R")
  source("step2.R")
  source("step3.R")
  source("Simulation 2/do_sim2.R")
  source("auxiliary_functions.R")
})

## load objects in cluster
export(backend, "cond_sim2")

## perform simulation
start  <- Sys.time()
output <- par_lapply(backend, 1:nrow(cond_sim2), do_sim2, cond = cond_sim2,
                     outputfile = "Simulation 2/sim2.csv", verbose = FALSE)
end <- Sys.time()

# close cluster:
stop_backend(backend)