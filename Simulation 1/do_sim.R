#### This script defines the simulation function do_sim() ####

do_sim <- function(pos, cond, outputfile, verbose = FALSE){
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished
  
  #### for testing:
  #pos = 1

  
  replication <- cond$replication[pos]
  iteration <- cond$iteration[pos]
  # get condition levels and set seed:
  n <- cond$n[pos]
  obs <- cond$obs[pos]
  invariance_level =  cond$invariance_level[pos] |> as.character()
  direction =  cond$direction[pos] |> as.character()
  seed_cond <- cond$seed[pos]
  set.seed(seed_cond)
  
  #### set data generation parameters ####
  ## regression parameters:
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
  
  ## innovation variances
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
  
  ## lambda matrix (loadings):
  loadings_baseline <- c(.55, .55, .7, .4)
  
  # loadings in group 1 are all equal to baseline:
  lambda_g1 <- lavaan::lav_matrix_bdiag(list(loadings_baseline,
                                             loadings_baseline))
  
  # loadings in group 2 are equal to group 1 unless modified by condition:
  loadings_g2 <- loadings_baseline
  if(invariance_level == "partial_metric"){
    if(direction == "mixed"){
      loadings_g2[3:4] <-  loadings_g2[4:3]
    } else {
      loadings_g2[3:4] <- loadings_g2[3:4] + .3
    }
  }
  lambda_g2 <- lavaan::lav_matrix_bdiag(list(loadings_g2,
                                             loadings_g2))
  
  ## theta matrix (residual variances):
  resvar_baseline <- c(.5, .9, .51, .84)
  
  # residual variances in group 1 are all equal to baseline:
  theta_g1 <- matrix(0, nrow = 8, ncol = 8)
  diag(theta_g1) <- c(resvar_baseline, resvar_baseline)
  
  # residual variances in group 2 are equal to group 1 unless modified by condition:
  resvar_g2 <- resvar_baseline
  if(invariance_level %in% c("full_scalar", "partial_scalar", "partial_metric")){
    if(direction == "mixed"){
      resvar_g2 <- resvar_g2[c(2, 1, 4, 3)]
    } else {
      resvar_g2 <- resvar_g2 - .4
    }
  }
  theta_g2 <- matrix(0, nrow = 8, ncol = 8)
  diag(theta_g2) <- c(resvar_g2, resvar_g2)
  
  ##  intercepts
  # intercepts are all 0 in group 1
  tau_g1 <- c(1, 1, 0, 2, 1, 1, 0, 2)
  
  # intercepts in group 2 are equal to group 1 unless modified by condition:
  tau_g2 <- tau_g1
  
  if(invariance_level %in% c("partial_scalar", "partial_metric")){
    if(direction == "mixed"){
      tau_g2[c(3, 4, 7, 8)] <- tau_g2[c(4, 3, 8, 7)]
    } else {
      tau_g2[c(3, 4, 7, 8)] <- tau_g2[c(3, 4, 7, 8)] + 1
    }
  }
  
  #### generate items scores ####
  # create empty data frame:
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
  
  ## create group assignment vector
  groupassignment <- c(rep("group1", n*0.5), rep("group2", n*0.5))
  
  for(i in 1:n){
    # create temporary dataframe:
    g <- groupassignment[i]
    
    # get correct SM and MM parameters:
    if(g == "group1"){
      phimat_i <- phimat_g1
      zetamat_i <- zetamat_g1
      lambda_i <- lambda_g1
      theta_i <- theta_g1
      tau_i <- tau_g1
      mu_i <- c(rnorm(1, 3.09, 1),
                rnorm(1, .98, 1))
    }
    if(g == "group2"){
      phimat_i <- phimat_g2
      zetamat_i <- zetamat_g2
      lambda_i <- lambda_g2
      theta_i <- theta_g2
      tau_i <- tau_g2
      mu_i <- c(rnorm(1, 4.56, 1),
                rnorm(1, .32, 1))
    }
    
    # generate factor scores:
    eta_i <- sim_VAR(factors = 2, obs = obs,
                     phi = phimat_i, zeta = zetamat_i,
                     mu = mu_i,
                     burn_in = 20)
    
    # generate errors:
    epsilon_i <- mvrnorm(obs, mu = rep(0, 8),
                         Sigma = theta_i, empirical=T)
    
    tau_matrix <- rep(tau_i, each = obs) |> matrix(nrow = obs)
    
    # transform factor scores into observed scores:
    data_i <- as.matrix(eta_i[, c("eta1", "eta2")]) %*% t(lambda_i) + tau_matrix + epsilon_i |>
      as.data.frame()
    colnames(data_i) <- paste0("v", 1:8)
    
    # add id and true cluster variable:
    data_i$id <- i
    data_i$group <- g
    data_i$obs <- eta_i$obs
    
    # merge with full data
    data <- dplyr::full_join(data, data_i, by = join_by(id, obs, group, v1, v2, v3, v4, v5, v6, v7, v8))
  }
  
  start <- Sys.time()
  model_step1 <- list(
    "f1 =~ 0.55*v1 + v2 + v3 + v4
      v1 ~ 1*1
      f1 ~ NA*1",
    "f2 =~ 0.55*v5 + v6 + v7 + v8
      v5 ~ 1*1
      f2 ~ NA*1")
  
  #### SINGLE GROUP ####
  #### Step 1 ####
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
                                  partial_noninvariances = partial_noninvariances,
                                  bounds = "wide")
  
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
    
    lambda_ests <- output_step2_multi$result$result$MMparameters$lambda_group
    lambda_ests <- do.call(rbind, lapply(lambda_ests, function(x){c(x[3:4, 1], x[7:8, 2])}))
    lambda_true <- matrix(c(lambda_g1[3:4, 1], lambda_g1[7:8, 2], 
                            lambda_g2[3:4, 1], lambda_g2[7:8, 2]),
                          nrow = 2, byrow = TRUE)
    bias_lambda <- sum(lambda_true - lambda_ests)/(4*2)  |> as.numeric()
    RMSE_lambda <- sqrt(sum((lambda_true - lambda_ests)^2)/(4*2)) |> as.numeric()
    
    theta_ests <- output_step2_multi$result$result$MMparameters$theta_group
    theta_ests <- do.call(rbind, lapply(theta_ests, diag))
    theta_ests <- theta_ests[, c(3, 4, 7, 8)]
    theta_true <- matrix(c(diag(theta_g1)[c(3, 4, 7, 8)],
                           diag(theta_g2)[c(3, 4, 7, 8)]),
                         nrow = 2, byrow = TRUE)
    bias_theta <- sum(theta_true - theta_ests)/(4*2)  |> as.numeric()
    RMSE_theta <- sqrt(sum((theta_true - theta_ests)^2)/(4*2)) |> as.numeric()
    
    tau_ests <- output_step2_multi$result$result$MMparameters$tau_group
    tau_ests <- do.call(rbind, lapply(tau_ests, function(x){c(x[3:4, 1], x[7:8, 2])}))
    tau_true <- matrix(c(tau_g1[c(3, 4, 7, 8)], tau_g2[c(3, 4, 7, 8)]),
                      nrow = 2, byrow = TRUE)
    bias_tau <- sum(tau_true - tau_ests)/(4*2)  |> as.numeric()
    RMSE_tau <- sqrt(sum((tau_true - tau_ests)^2)/(4*2)) |> as.numeric()
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
    
    bias_lambda <- NA
    RMSE_lambda <- NA
    bias_theta <- NA
    RMSE_theta <- NA
    bias_tau <- NA
    RMSE_tau <- NA
  }
  
  duration <- difftime(Sys.time(), start, units = "secs") |> as.numeric()
  
  output <- c("iteration" = iteration, "replication" = replication,
              "n" = n, "obs" = obs,
              "invariance_level" = invariance_level,
              "direction" = direction,
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
              "bias_lambda" = bias_lambda, "bias_theta" = bias_theta, "bias_tau" = bias_tau,
              "RMSE_lambda" = RMSE_lambda, "RMSE_theta" = RMSE_theta, "RMSE_tau" = RMSE_tau,
              "step1_single_warning" = step1_single_warning, "step2_single_warning" = step2_single_warning, "step3_single_warning" = step3_single_warning,
              "step1_single_error" = step1_single_error, "step2_single_error" = step2_single_error, "step3_single_error" = step3_single_error,
              "step1_multi_warning" = step1_multi_warning, "step2_multi_warning" = step2_multi_warning, "step3_multi_warning" = step3_multi_warning,
              "step1_multi_error" = step1_multi_error, "step2_multi_error" = step2_multi_error, "step3_multi_error" = step3_multi_error,
              "seed" = seed_cond, "pos" = pos,
              "step1_single_warning_text" = step1_single_warning_text, "step2_single_warning_text" = step2_single_warning_text, "step3_single_warning_text" = step3_single_warning_text,
              "step1_single_error_text" = step1_single_error_text, "step2_single_error_text" = step2_single_error_text, "step3_single_error_text" = step3_single_error_text,
              "step1_multi_warning_text" = step1_multi_warning_text, "step2_multi_warning_text" = step2_multi_warning_text, "step3_multi_warning_text" = step3_multi_warning_text,
              "step1_multi_error_text" = step1_multi_error_text, "step2_multi_error_text" = step2_multi_error_text, "step3_multi_error_text" = step3_multi_error_text)
  
  for(i in 98:109){
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