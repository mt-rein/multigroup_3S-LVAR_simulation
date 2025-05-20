#### This script defines the simulation function do_sim() ####

do_sim <- function(pos, cond, outputfile, verbose = FALSE){
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished
  
  #### for testing:
  # pos = 1
  # invariance_level <- "partial_metric"
  # pattern = "unidirectional"
  # ss_n = 100
  # ss_t = 40
  # ss_ratio = "balanced"

  
  replication <- cond$replication[pos]
  iteration <- cond$iteration[pos]
  # get condition levels and set seed:
  invariance_level <- cond$invariance_level[pos] |> as.character()
  pattern <- cond$pattern[pos] |> as.character()
  ss_n <- cond$ss_n[pos]
  ss_t <- cond$ss_t[pos]
  ss_ratio = cond$ss_ratio[pos] |> as.character()
  seed_cond <- cond$seed[pos]
  set.seed(seed_cond)
  
  #### set data generation parameters ####
  ## structural parameters
  ## regression coefficients:
  # group 1
  phi11_g1_pop <- .33
  phi22_g1_pop <- .37
  phi12_g1_pop <- .05
  phi21_g1_pop <- .04
  
  # group 2
  phi11_g2_pop <- .42
  phi22_g2_pop <- .23
  phi12_g2_pop <- .17
  phi21_g2_pop <- .02
  
  # combine into matrices:
  phimat_g1 <- matrix(c(phi11_g1_pop, phi12_g1_pop,
                        phi21_g1_pop, phi22_g1_pop),
                      ncol = 2, byrow = TRUE)
  phimat_g2 <- matrix(c(phi11_g2_pop, phi12_g2_pop,
                        phi21_g2_pop, phi22_g2_pop),
                      ncol = 2, byrow = TRUE)
  
  ## innovation variances:
  # total covariances are equal in both groups:
  psimat <- matrix(c(1, .3, .3, 1), ncol = 2)
  
  # innovation covariances depend on group-specific regression effects:
  zetamat_g1 <- (diag(2*2) - kronecker(phimat_g1, phimat_g1)) %*% c(psimat) |> matrix(ncol = 2)
  zeta1_g1_pop <- zetamat_g1[1, 1]
  zeta2_g1_pop <- zetamat_g1[2, 2]
  zeta12_g1_pop <- zetamat_g1[1, 2]
  zetamat_g2 <- (diag(2*2) - kronecker(phimat_g2, phimat_g2)) %*% c(psimat) |> matrix(ncol = 2)
  zeta1_g2_pop <- zetamat_g2[1, 1]
  zeta2_g2_pop <- zetamat_g2[2, 2]
  zeta12_g2_pop <- zetamat_g2[1, 2]
  
  
  ## measurement parameters:
  # baseline:
  loadings_g1 <- loadings_g2 <- c(.75, .75, .5, .5)
  resvar_g1 <- 1-(loadings_g1^2)
  resvar_g2 <- 1-(loadings_g2^2)
  tau_g1 <- tau_g2 <- c(1, 1, 1, 1)
  
  ## depending on condition, introduce differences in residual variances, intercepts, and loadings
  # residual variances:
  if(invariance_level %in% c("full_scalar", "partial_scalar", "partial_metric")){
    if(pattern == "mixed"){
      resvar_g1[c(1, 3)] <- resvar_g1[c(1, 3)] - .25
      resvar_g1[c(2, 4)] <- resvar_g1[c(2, 4)] + .25
      resvar_g2[c(1, 3)] <- resvar_g2[c(1, 3)] + .25
      resvar_g2[c(2, 4)] <- resvar_g2[c(2, 4)] - .25
    }
    if(pattern == "unidirectional"){
      resvar_g1 <- resvar_g1 - .25
      resvar_g2 <- resvar_g2 + .25
    }
  }
  
  theta_g1 <- diag(c(resvar_g1, resvar_g1))
  theta_g2 <- diag(c(resvar_g2, resvar_g2))
  
  # intercepts:
  if(invariance_level %in% c("partial_scalar", "partial_metric")){
    if(pattern == "mixed"){
      tau_g1[3] <-  tau_g1[3] + .5
      tau_g1[4] <-  tau_g1[4] - .5
      tau_g2[3] <-  tau_g2[3] - .5
      tau_g2[4] <-  tau_g2[4] + .5
    }
    if(pattern == "unidirectional"){
      tau_g1[3:4] <- tau_g1[3:4] + .5
      tau_g2[3:4] <- tau_g2[3:4] - .5
    }
  }
  tau_g1 <- c(tau_g1, tau_g1)
  tau_g2 <- c(tau_g2, tau_g2)
  
  # loadings:
  if(invariance_level == "partial_metric"){
    if(pattern == "mixed"){
      loadings_g1[3] <- loadings_g1[3] + .15
      loadings_g1[4] <- loadings_g1[4] - .15
      loadings_g2[3] <- loadings_g2[3] - .15
      loadings_g2[4] <- loadings_g2[4] + .15
    }
    if(pattern == "unidirectional"){
      loadings_g1[3:4] <- loadings_g1[3:4] + .15
      loadings_g2[3:4] <- loadings_g2[3:4] - .15
    }
  }
  lambda_g1 <- lavaan::lav_matrix_bdiag(list(loadings_g1, loadings_g1))
  lambda_g2 <- lavaan::lav_matrix_bdiag(list(loadings_g2, loadings_g2))
  
  
  ## calculate population reliability:
  sigma_g1 <- lambda_g1 %*% psimat %*% t(lambda_g1) + theta_g1
  sigma_g2 <- lambda_g2 %*% psimat %*% t(lambda_g2) + theta_g2
  lambdastar_g1_pop <- diag(psimat %*% t(lambda_g1) %*% solve(sigma_g1) %*% lambda_g1)
  lambdastar_g2_pop <- diag(psimat %*% t(lambda_g2) %*% solve(sigma_g2) %*% lambda_g2)
  lambdastar_f1_g1_pop <- lambdastar_g1_pop[1]
  lambdastar_f2_g1_pop <- lambdastar_g1_pop[2]
  lambdastar_f1_g2_pop <- lambdastar_g2_pop[1]
  lambdastar_f2_g2_pop <- lambdastar_g2_pop[2]
  
  #### generate observed items scores ####
  # sample size per group
  if(ss_ratio == "balanced"){
    n1 <- n2 <- ss_n*0.5
  }
  if(ss_ratio == "unbalanced"){
    n1 <- ss_n*0.25
    n2 <- ss_n*0.75
  }
  
  # create empty data frame for all (observed) items
  data <- data.frame(id = numeric(),
                     obs = numeric(),
                     group = character(),
                     v1 = numeric(),
                     v2 = numeric(),
                     v3 = numeric(),
                     v4 = numeric(),
                     v5 = numeric(),
                     v6 = numeric(),
                     v7 = numeric(),
                     v8 = numeric())
  
  ## generate factor scores and the observed items per group:
  for(g in c("group1", "group2")){
    # create empty data frame for group's factor scores:
    eta_g <- data.frame(id = numeric(),
                        obs = numeric(),
                        eta1 = numeric(),
                        eta2 = numeric())
    
    # get correct SM and MM parameters:
    if(g == "group1"){
      ids_g <- 1:n1
      phimat_g <- phimat_g1
      zetamat_g <- zetamat_g1
      mu_g <- c(3.09, .98)
      lambda_g <- lambda_g1
      theta_g <- theta_g1
      tau_g <- tau_g1
    } else {
      ids_g <- (1:n2)+n1
      phimat_g <- phimat_g2
      zetamat_g <- zetamat_g2
      mu_g <- c(4.56, .32)
      lambda_g <- lambda_g2
      theta_g <- theta_g2
      tau_g <- tau_g2
    }
    for(i in ids_g){
      # generate person-specific latent means:
      mu_i <- c(rnorm(1, mu_g[1], 1),
                rnorm(1, mu_g[2], 1))
      
      # generate time series data (factor scores) for individual:
      eta_i <- sim_VAR(factors = 2, obs = ss_t,
                       phi = phimat_g, zeta = zetamat_g,
                       mu = mu_i,
                       burn_in = 20)
      
      eta_i$id <- i
      # merge with group's factor score data frame:
      eta_g <- dplyr::full_join(eta_g, eta_i, by = join_by(id, obs, eta1, eta2))
    }
    
    # generate errors (item residuals):
    epsilon_g <- mvrnorm(nrow(eta_g), mu = rep(0, 8),
                       Sigma = theta_g, empirical = FALSE)
    # transform factor scores into group's observed scores:
    data_g <- t(tau_g + lambda_g %*% t(eta_g[, c("eta1", "eta2")])) + epsilon_g |>
      as.data.frame()
    colnames(data_g) <- paste0("v", 1:8)
    data_g$id <- eta_g$id
    data_g$obs <- eta_g$obs
    data_g$group <- g
    
    # merge group's observed items with full data frame:
    data <- dplyr::full_join(data, data_g, by = join_by(id, obs, group, v1, v2, v3, v4, v5, v6, v7, v8))
  }
  
  start <- Sys.time()
  model_step1 <- list(
    "f1 =~ 0.75*v1 + v2 + v3 + v4
      v1 ~ 1*1
      f1 ~ NA*1",
    "f2 =~ 0.75*v5 + v6 + v7 + v8
      v5 ~ 1*1
      f2 ~ NA*1")
  
  #### SINGLE GROUP ####
  #### Step 1 ####
  output_step1_single <- run_step1(data = data,
                                   measurementmodel = model_step1,
                                   group = NULL)
  
  # extract error/warning messages (if applicable):
  step1_single_warning <- ifelse(is_empty(output_step1_single$warnings),
                                 FALSE, TRUE)
  step1_single_warning_text <- ifelse(is_empty(output_step1_single$warnings),
                                      "",
                                      paste(c(output_step1_single$warnings),
                                            collapse = "; ")
  )
  step1_single_error <- ifelse(is_empty(output_step1_single$result$error),
                               FALSE, TRUE)
  step1_single_error_text <- ifelse(is_empty(output_step1_single$result$error),
                                    "",
                                    paste(c(output_step1_single$result$error),
                                          collapse = "; "))
  
  # check whether the solution is admissible:
  if(grepl("variances are negative", step1_single_warning_text, ignore.case = TRUE)) {
    # Re-run step 1 with wide bounds in case there's a heywood case
    output_step1_single <- run_step1(data = data,
                                     measurementmodel = model_step1,
                                     group = NULL,
                                     bounds = "wide")
    
    # extract error/warning messages (if applicable):
    step1_single_warning <- ifelse(is_empty(output_step1_single$warnings),
                                   FALSE, TRUE)
    step1_single_warning_text <- ifelse(is_empty(output_step1_single$warnings),
                                        "",
                                        paste(c(output_step1_single$warnings),
                                              collapse = "; ")
    )
    step1_single_error <- ifelse(is_empty(output_step1_single$result$error),
                                 FALSE, TRUE)
    step1_single_error_text <- ifelse(is_empty(output_step1_single$result$error),
                                      "",
                                      paste(c(output_step1_single$result$error),
                                            collapse = "; "))
    rerun_step1 <- TRUE
  } else {
    rerun_step1 <- FALSE
  }
  
  #### Step 2 ####
  if(!step1_single_error){                                                             # only proceed if there is no error in step 1
    output_step2_single <- run_step2(step1output = output_step1_single$result$result)
    # extract error/warning messages (if applicable):
    step2_single_warning <- ifelse(is_empty(output_step2_single$warnings),
                                   FALSE, TRUE)
    step2_single_warning_text <- ifelse(is_empty(output_step2_single$warnings),
                                        "",
                                        paste(c(output_step2_single$warnings),
                                              collapse = "; ")
    )
    step2_single_error <- ifelse(is_empty(output_step2_single$result$error),
                                 FALSE, TRUE)
    step2_single_error_text <- ifelse(is_empty(output_step2_single$result$error),
                                      "",
                                      paste(c(output_step2_single$result$error),
                                            collapse = "; ")
    )
  } else {
    step2_single_warning <- FALSE
    step2_single_warning_text <- "step1 not successful"
    step2_single_error <- FALSE
    step2_single_error_text <- "step1 not successful"
  }
  
  #### Step 3 ####
  if(!step1_single_error & !step2_single_error){                                              # only proceed if there is no error in step 1 as well as step 2
    output_step3_single <- run_step3(step2output = output_step2_single$result$result,
                                     id = "id", step3group = "group")
    
    # extract error/warning messages (if applicable):
    step3_single_warning <- ifelse(is_empty(output_step3_single$warnings),
                                   FALSE, TRUE)
    step3_single_warning_text <- ifelse(is_empty(output_step3_single$warnings),
                                        "",
                                        paste(c(output_step3_single$warnings),
                                              collapse = "; ")
    )
    step3_single_error <- ifelse(is_empty(output_step3_single$result$error),
                                 FALSE, TRUE)
    step3_single_error_text <- ifelse(is_empty(output_step3_single$result$error),
                                      "",
                                      paste(c(output_step3_single$result$error),
                                            collapse = "; ")
    )
    
    # check if the model converged:
    if(output_step3_single$result$result$model@output$status$code != 0){
      step3_single_error <- TRUE
      step3_single_error_text <- "step3 model estimation failed"
    }
  } else {
    step3_single_warning <- FALSE
    step3_single_warning_text <- "step1 or step2 not successful"
    step3_single_error <- FALSE
    step3_single_error_text <- "step1 or step2 not successful"
  }
  
  if(!step1_single_error & !step2_single_error & !step3_single_error){
    estimates <- output_step3_single$result$result$estimates
    standarderrors <- setNames(summary(output_step3_single$result$result$model)$parameters$Std.Error,
                               summary(output_step3_single$result$result$model)$parameters$name)
    
    single_phi11_g1 <- estimates["phi11_group1"] |> as.numeric()
    single_phi22_g1 <- estimates["phi22_group1"] |> as.numeric()
    single_phi12_g1 <- estimates["phi12_group1"] |> as.numeric()
    single_phi21_g1 <- estimates["phi21_group1"] |> as.numeric()
    
    single_zeta1_g1 <- estimates["zeta1_group1"] |> as.numeric()
    single_zeta2_g1 <- estimates["zeta2_group1"] |> as.numeric()
    single_zeta12_g1 <- estimates["zeta12_group1"] |> as.numeric()
    
    single_phi11_g2 <- estimates["phi11_group2"] |> as.numeric()
    single_phi22_g2 <- estimates["phi22_group2"] |> as.numeric()
    single_phi12_g2 <- estimates["phi12_group2"] |> as.numeric()
    single_phi21_g2 <- estimates["phi21_group2"] |> as.numeric()
    
    single_zeta1_g2 <- estimates["zeta1_group2"] |> as.numeric()
    single_zeta2_g2 <- estimates["zeta2_group2"] |> as.numeric()
    single_zeta12_g2 <- estimates["zeta12_group2"] |> as.numeric()
    
    single_phi11_g1_se <- standarderrors["phi11_group1"] |> as.numeric()
    single_phi22_g1_se <- standarderrors["phi22_group1"] |> as.numeric()
    single_phi12_g1_se <- standarderrors["phi12_group1"] |> as.numeric()
    single_phi21_g1_se <- standarderrors["phi21_group1"] |> as.numeric()
    
    single_zeta1_g1_se <- standarderrors["zeta1_group1"] |> as.numeric()
    single_zeta2_g1_se <- standarderrors["zeta2_group1"] |> as.numeric()
    single_zeta12_g1_se <- standarderrors["zeta12_group1"] |> as.numeric()
    
    single_phi11_g2_se <- standarderrors["phi11_group2"] |> as.numeric()
    single_phi22_g2_se <- standarderrors["phi22_group2"] |> as.numeric()
    single_phi12_g2_se <- standarderrors["phi12_group2"] |> as.numeric()
    single_phi21_g2_se <- standarderrors["phi21_group2"] |> as.numeric()
    
    single_zeta1_g2_se <- standarderrors["zeta1_group2"] |> as.numeric()
    single_zeta2_g2_se <- standarderrors["zeta2_group2"] |> as.numeric()
    single_zeta12_g2_se <- standarderrors["zeta12_group2"] |> as.numeric()
    
  } else {
    single_phi11_g1 <- NA
    single_phi22_g1 <- NA
    single_phi12_g1 <- NA
    single_phi21_g1 <- NA
    
    single_zeta1_g1 <- NA
    single_zeta2_g1 <- NA
    single_zeta12_g1 <- NA
    
    single_phi11_g2 <- NA
    single_phi22_g2 <- NA
    single_phi12_g2 <- NA
    single_phi21_g2 <- NA
    
    single_zeta1_g2 <- NA
    single_zeta2_g2 <- NA
    single_zeta12_g2 <- NA
    
    single_phi11_g1_se <- NA
    single_phi22_g1_se <- NA
    single_phi12_g1_se <- NA
    single_phi21_g1_se <- NA
    
    single_zeta1_g1_se <- NA
    single_zeta2_g1_se <- NA
    single_zeta12_g1_se <- NA
    
    single_phi11_g2_se <- NA
    single_phi22_g2_se <- NA
    single_phi12_g2_se <- NA
    single_phi21_g2_se <- NA
    
    single_zeta1_g2_se <- NA
    single_zeta2_g2_se <- NA
    single_zeta12_g2_se <- NA
  }
  
  
  #### MULTI-GROUP ####
  # which parameters are constrained to equality?
  if(invariance_level == "full_strict"){
    invariances <- c("loadings", "intercepts", "residuals")
  } else {
    invariances <- c("loadings", "intercepts")
  }
  
  partial_noninvariances <- list(NULL, NULL)
  # if there are partial invariances, add these:
  if(invariance_level == "partial_scalar"){
    partial_noninvariances <- list(c("v3 ~ 1", "v4 ~ 1"),
                                   c("v7 ~ 1", "v8 ~ 1"))
  }
  if(invariance_level == "partial_metric"){
    partial_noninvariances <- list(c("f1 =~ v3", "f1 =~ v4",
                                     "v3 ~ 1", "v4 ~ 1"),
                                   c("f2 =~ v7", "f2 =~ v8",
                                     "v7 ~ 1", "v8 ~ 1"))
  }
  
  output_step1_multi <- run_step1(data = data,
                                  measurementmodel = model_step1,
                                  group = "group",
                                  invariances = invariances,
                                  partial_noninvariances = partial_noninvariances)
  
  # extract error/warning messages (if applicable):
  step1_multi_warning <- ifelse(is_empty(output_step1_multi$warnings),
                                FALSE, TRUE)
  step1_multi_warning_text <- ifelse(is_empty(output_step1_multi$warnings),
                                     "",
                                     paste(c(output_step1_multi$warnings),
                                           collapse = "; ")
  )
  step1_multi_error <- ifelse(is_empty(output_step1_multi$result$error),
                              FALSE, TRUE)
  step1_multi_error_text <- ifelse(is_empty(output_step1_multi$result$error),
                                   "",
                                   paste(c(output_step1_multi$result$error),
                                         collapse = "; "))
  
  #### Step 2 ####
  if(!step1_multi_error){                                                             # only proceed if there is no error in step 1
    output_step2_multi <- run_step2(step1output = output_step1_multi$result$result)
    # extract error/warning messages (if applicable):
    step2_multi_warning <- ifelse(is_empty(output_step2_multi$warnings),
                                  FALSE, TRUE)
    step2_multi_warning_text <- ifelse(is_empty(output_step2_multi$warnings),
                                       "",
                                       paste(c(output_step2_multi$warnings),
                                             collapse = "; ")
    )
    step2_multi_error <- ifelse(is_empty(output_step2_multi$result$error),
                                FALSE, TRUE)
    step2_multi_error_text <- ifelse(is_empty(output_step2_multi$result$error),
                                     "",
                                     paste(c(output_step2_multi$result$error),
                                           collapse = "; ")
    )
  } else {
    step2_multi_warning <- FALSE
    step2_multi_warning_text <- "step1 not successful"
    step2_multi_error <- FALSE
    step2_multi_error_text <- "step1 not successful"
  }
  
  #### Step 3 ####
  if(!step1_multi_error & !step2_multi_error){                                              # only proceed if there is no error in step 1 as well as step 2
    output_step3_multi <- run_step3(step2output = output_step2_multi$result$result,
                                    id = "id", step3group = "group")
    
    # extract error/warning messages (if applicable):
    step3_multi_warning <- ifelse(is_empty(output_step3_multi$warnings),
                                  FALSE, TRUE)
    step3_multi_warning_text <- ifelse(is_empty(output_step3_multi$warnings),
                                       "",
                                       paste(c(output_step3_multi$warnings),
                                             collapse = "; ")
    )
    step3_multi_error <- ifelse(is_empty(output_step3_multi$result$error),
                                FALSE, TRUE)
    step3_multi_error_text <- ifelse(is_empty(output_step3_multi$result$error),
                                     "",
                                     paste(c(output_step3_multi$result$error),
                                           collapse = "; ")
    )
    
    # check if the model converged:
    if(output_step3_multi$result$result$model@output$status$code != 0){
      step3_multi_error <- TRUE
      step3_multi_error_text <- "step3 model estimation failed"
    }
    
    
  } else {
    step3_multi_warning <- FALSE
    step3_multi_warning_text <- "step1 or step2 not successful"
    step3_multi_error <- FALSE
    step3_multi_error_text <- "step1 or step2 not successful"
  }
  
  if(!step1_multi_error & !step2_multi_error & !step3_multi_error){
    estimates <- output_step3_multi$result$result$estimates
    standarderrors <- setNames(summary(output_step3_multi$result$result$model)$parameters$Std.Error,
                               summary(output_step3_multi$result$result$model)$parameters$name)
    
    multi_phi11_g1 <- estimates["phi11_group1"] |> as.numeric()
    multi_phi22_g1 <- estimates["phi22_group1"] |> as.numeric()
    multi_phi12_g1 <- estimates["phi12_group1"] |> as.numeric()
    multi_phi21_g1 <- estimates["phi21_group1"] |> as.numeric()
    
    multi_zeta1_g1 <- estimates["zeta1_group1"] |> as.numeric()
    multi_zeta2_g1 <- estimates["zeta2_group1"] |> as.numeric()
    multi_zeta12_g1 <- estimates["zeta12_group1"] |> as.numeric()
    
    multi_phi11_g2 <- estimates["phi11_group2"] |> as.numeric()
    multi_phi22_g2 <- estimates["phi22_group2"] |> as.numeric()
    multi_phi12_g2 <- estimates["phi12_group2"] |> as.numeric()
    multi_phi21_g2 <- estimates["phi21_group2"] |> as.numeric()
    
    multi_zeta1_g2 <- estimates["zeta1_group2"] |> as.numeric()
    multi_zeta2_g2 <- estimates["zeta2_group2"] |> as.numeric()
    multi_zeta12_g2 <- estimates["zeta12_group2"] |> as.numeric()
    
    multi_phi11_g1_se <- standarderrors["phi11_group1"] |> as.numeric()
    multi_phi22_g1_se <- standarderrors["phi22_group1"] |> as.numeric()
    multi_phi12_g1_se <- standarderrors["phi12_group1"] |> as.numeric()
    multi_phi21_g1_se <- standarderrors["phi21_group1"] |> as.numeric()
    
    multi_zeta1_g1_se <- standarderrors["zeta1_group1"] |> as.numeric()
    multi_zeta2_g1_se <- standarderrors["zeta2_group1"] |> as.numeric()
    multi_zeta12_g1_se <- standarderrors["zeta12_group1"] |> as.numeric()
    
    multi_phi11_g2_se <- standarderrors["phi11_group2"] |> as.numeric()
    multi_phi22_g2_se <- standarderrors["phi22_group2"] |> as.numeric()
    multi_phi12_g2_se <- standarderrors["phi12_group2"] |> as.numeric()
    multi_phi21_g2_se <- standarderrors["phi21_group2"] |> as.numeric()
    
    multi_zeta1_g2_se <- standarderrors["zeta1_group2"] |> as.numeric()
    multi_zeta2_g2_se <- standarderrors["zeta2_group2"] |> as.numeric()
    multi_zeta12_g2_se <- standarderrors["zeta12_group2"] |> as.numeric()
    
    lambdastar_f1_g1 <- output_step2_multi$result$result$lambda_star[1, 1]
    lambdastar_f2_g1 <- output_step2_multi$result$result$lambda_star[1, 2]
    lambdastar_f1_g2 <- output_step2_multi$result$result$lambda_star[2, 1]
    lambdastar_f2_g2 <- output_step2_multi$result$result$lambda_star[2, 2]
    
    lambda_ests <- output_step2_multi$result$result$MMparameters$lambda_group
    lambda_ests <- do.call(rbind, lapply(lambda_ests, function(x){c(x[3:4, 1], x[7:8, 2])}))
    lambda_true <- matrix(c(lambda_g1[3:4, 1], lambda_g1[7:8, 2], 
                            lambda_g2[3:4, 1], lambda_g2[7:8, 2]),
                          nrow = 2, byrow = TRUE)
    bias_lambda <- sum(lambda_ests - lambda_true)/(4*2)  |> as.numeric()
    RMSE_lambda <- sqrt(sum((lambda_ests - lambda_true)^2)/(4*2)) |> as.numeric()
    
    theta_ests <- output_step2_multi$result$result$MMparameters$theta_group
    theta_ests <- do.call(rbind, lapply(theta_ests, diag))
    theta_true <- matrix(c(diag(theta_g1), diag(theta_g2)),
                         nrow = 2, byrow = TRUE)
    bias_theta <- sum(theta_ests - theta_true)/(8*2)  |> as.numeric()
    RMSE_theta <- sqrt(sum((theta_ests - theta_true)^2)/(8*2)) |> as.numeric()
    
    tau_ests <- output_step2_multi$result$result$MMparameters$tau_group
    tau_ests <- do.call(rbind, lapply(tau_ests, function(x){c(x[3:4, 1], x[7:8, 2])}))
    tau_true <- matrix(c(tau_g1[c(3, 4, 7, 8)], tau_g2[c(3, 4, 7, 8)]),
                      nrow = 2, byrow = TRUE)
    bias_tau <- sum(tau_ests - tau_true)/(4*2)  |> as.numeric()
    RMSE_tau <- sqrt(sum((tau_ests - tau_true)^2)/(4*2)) |> as.numeric()
  } else {
    multi_phi11_g1 <- NA
    multi_phi22_g1 <- NA
    multi_phi12_g1 <- NA
    multi_phi21_g1 <- NA
    
    multi_zeta1_g1 <- NA
    multi_zeta2_g1 <- NA
    multi_zeta12_g1 <- NA
    
    multi_phi11_g2 <- NA
    multi_phi22_g2 <- NA
    multi_phi12_g2 <- NA
    multi_phi21_g2 <- NA
    
    multi_zeta1_g2 <- NA
    multi_zeta2_g2 <- NA
    multi_zeta12_g2 <- NA
    
    multi_phi11_g1_se <- NA
    multi_phi22_g1_se <- NA
    multi_phi12_g1_se <- NA
    multi_phi21_g1_se <- NA
    
    multi_zeta1_g1_se <- NA
    multi_zeta2_g1_se <- NA
    multi_zeta12_g1_se <- NA
    
    multi_phi11_g2_se <- NA
    multi_phi22_g2_se <- NA
    multi_phi12_g2_se <- NA
    multi_phi21_g2_se <- NA
    
    multi_zeta1_g2_se <- NA
    multi_zeta2_g2_se <- NA
    multi_zeta12_g2_se <- NA
    
    lambdastar_f1_g1 <- NA
    lambdastar_f2_g1 <- NA
    lambdastar_f1_g2 <- NA
    lambdastar_f2_g2 <- NA
    
    bias_lambda <- NA
    RMSE_lambda <- NA
    bias_theta <- NA
    RMSE_theta <- NA
    bias_tau <- NA
    RMSE_tau <- NA
  }
  
  duration <- difftime(Sys.time(), start, units = "secs") |> as.numeric()
  
  output <- c("iteration" = iteration, "replication" = replication,
              "invariance_level" = invariance_level, "pattern" = pattern,
              "ss_n" = ss_n, "ss_t" = ss_t, "ss_ratio" = ss_ratio,
              "duration" = duration,
              "phi11_g1_pop" = phi11_g1_pop, "phi12_g1_pop" = phi12_g1_pop, "phi21_g1_pop" = phi21_g1_pop, "phi22_g1_pop" = phi22_g1_pop,
              "phi11_g2_pop" = phi11_g2_pop, "phi12_g2_pop" = phi12_g2_pop, "phi21_g2_pop" = phi21_g2_pop, "phi22_g2_pop" = phi22_g2_pop,
              "zeta1_g1_pop" = zeta1_g1_pop, "zeta2_g1_pop" = zeta2_g1_pop, "zeta12_g1_pop" = zeta12_g1_pop,
              "zeta1_g2_pop" = zeta1_g2_pop, "zeta2_g2_pop" = zeta2_g2_pop, "zeta12_g2_pop" = zeta12_g2_pop,
              "single_phi11_g1" = single_phi11_g1, "single_phi12_g1" = single_phi12_g1, "single_phi21_g1" = single_phi21_g1, "single_phi22_g1" = single_phi22_g1,
              "single_phi11_g2" = single_phi11_g2, "single_phi12_g2" = single_phi12_g2, "single_phi21_g2" = single_phi21_g2, "single_phi22_g2" = single_phi22_g2,
              "single_zeta1_g1" = single_zeta1_g1, "single_zeta2_g1" = single_zeta2_g1, "single_zeta12_g1" = single_zeta12_g1,
              "single_zeta1_g2" = single_zeta1_g2, "single_zeta2_g2" = single_zeta2_g2, "single_zeta12_g2" = single_zeta12_g2,
              "single_phi11_g1_se" = single_phi11_g1_se, "single_phi12_g1_se" = single_phi12_g1_se, "single_phi21_g1_se" = single_phi21_g1_se, "single_phi22_g1_se" = single_phi22_g1_se,
              "single_phi11_g2_se" = single_phi11_g2_se, "single_phi12_g2_se" = single_phi12_g2_se, "single_phi21_g2_se" = single_phi21_g2_se, "single_phi22_g2_se" = single_phi22_g2_se,
              "single_zeta1_g1_se" = single_zeta1_g1_se, "single_zeta2_g1_se" = single_zeta2_g1_se, "single_zeta12_g1_se" = single_zeta12_g1_se,
              "single_zeta1_g2_se" = single_zeta1_g2_se, "single_zeta2_g2_se" = single_zeta2_g2_se, "single_zeta12_g2_se" = single_zeta12_g2_se,
              "multi_phi11_g1" = multi_phi11_g1, "multi_phi12_g1" = multi_phi12_g1, "multi_phi21_g1" = multi_phi21_g1, "multi_phi22_g1" = multi_phi22_g1,
              "multi_phi11_g2" = multi_phi11_g2, "multi_phi12_g2" = multi_phi12_g2, "multi_phi21_g2" = multi_phi21_g2, "multi_phi22_g2" = multi_phi22_g2,
              "multi_zeta1_g1" = multi_zeta1_g1, "multi_zeta2_g1" = multi_zeta2_g1, "multi_zeta12_g1" = multi_zeta12_g1,
              "multi_zeta1_g2" = multi_zeta1_g2, "multi_zeta2_g2" = multi_zeta2_g2, "multi_zeta12_g2" = multi_zeta12_g2,
              "multi_phi11_g1_se" = multi_phi11_g1_se, "multi_phi12_g1_se" = multi_phi12_g1_se, "multi_phi21_g1_se" = multi_phi21_g1_se, "multi_phi22_g1_se" = multi_phi22_g1_se,
              "multi_phi11_g2_se" = multi_phi11_g2_se, "multi_phi12_g2_se" = multi_phi12_g2_se, "multi_phi21_g2_se" = multi_phi21_g2_se, "multi_phi22_g2_se" = multi_phi22_g2_se,
              "multi_zeta1_g1_se" = multi_zeta1_g1_se, "multi_zeta2_g1_se" = multi_zeta2_g1_se, "multi_zeta12_g1_se" = multi_zeta12_g1_se,
              "multi_zeta1_g2_se" = multi_zeta1_g2_se, "multi_zeta2_g2_se" = multi_zeta2_g2_se, "multi_zeta12_g2_se" = multi_zeta12_g2_se,
              "lambdastar_f1_g1_pop" = lambdastar_f1_g1_pop, "lambdastar_f2_g1_pop" = lambdastar_f2_g1_pop,
              "lambdastar_f1_g2_pop" = lambdastar_f1_g2_pop, "lambdastar_f2_g2_pop" = lambdastar_f2_g2_pop,
              "lambdastar_f1_g1" = lambdastar_f1_g1, "lambdastar_f2_g1" =  lambdastar_f2_g1,
              "lambdastar_f1_g2" = lambdastar_f1_g2, "lambdastar_f2_g2" =  lambdastar_f2_g2,
              "bias_lambda" = bias_lambda, "bias_theta" = bias_theta, "bias_tau" = bias_tau,
              "RMSE_lambda" = RMSE_lambda, "RMSE_theta" = RMSE_theta, "RMSE_tau" = RMSE_tau,
              "step1_single_warning" = step1_single_warning, "step2_single_warning" = step2_single_warning, "step3_single_warning" = step3_single_warning,
              "step1_single_error" = step1_single_error, "step2_single_error" = step2_single_error, "step3_single_error" = step3_single_error,
              "rerun_step1" = rerun_step1,
              "step1_multi_warning" = step1_multi_warning, "step2_multi_warning" = step2_multi_warning, "step3_multi_warning" = step3_multi_warning,
              "step1_multi_error" = step1_multi_error, "step2_multi_error" = step2_multi_error, "step3_multi_error" = step3_multi_error,
              "seed" = seed_cond, "pos" = pos,
              "step1_single_warning_text" = step1_single_warning_text, "step2_single_warning_text" = step2_single_warning_text, "step3_single_warning_text" = step3_single_warning_text,
              "step1_single_error_text" = step1_single_error_text, "step2_single_error_text" = step2_single_error_text, "step3_single_error_text" = step3_single_error_text,
              "step1_multi_warning_text" = step1_multi_warning_text, "step2_multi_warning_text" = step2_multi_warning_text, "step3_multi_warning_text" = step3_multi_warning_text,
              "step1_multi_error_text" = step1_multi_error_text, "step2_multi_error_text" = step2_multi_error_text, "step3_multi_error_text" = step3_multi_error_text)
  
  for(i in 108:119){
    output[i] <- str_squish(output[i])                                          # removes all whitespace and linebreaks from the error and warning strings
    output[i] <- gsub(",", "", output[i])                                       # removes all commata from error and warning strings (to prevent messing up the CSV file)
  }
  
  # check if file exists
  if(!file.exists(outputfile)){
    # if file does not yet exist
    write.table(t(output), file = outputfile, append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    # lock the file to prevent multiple processes accessing it simultaneously
    lock <- flock::lock(outputfile)
    write.table(t(output), file = outputfile, append = TRUE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
    # unlock the file
    flock::unlock(lock)
  }
  
  if(verbose == TRUE){
    print(paste("Simulation", pos, "completed at", Sys.time()))                 # prints a message when a replication is done (as a sign that R did not crash)
  }
  
  return(output)
}