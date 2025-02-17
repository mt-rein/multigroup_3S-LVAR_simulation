#### This script defines the simulation function do_sim() ####

do_sim <- function(pos, cond, outputfile, verbose = FALSE){
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished
  
  #### for testing:
  # pos = 1
  # replication <- 1
  # iteration <- 1
  # # get condition levels and set seed:
  # n <- 50
  # obs <- 25
  # lambda_noninvariance =  "uniform"
  # theta_noninvariance =  "uniform"
  # nu_noninvariance =  "uniform"
  # seed_cond <- 12345
  # set.seed(seed_cond)
  
  replication <- cond$replication[pos]
  iteration <- cond$iteration[pos]
  # get condition levels and set seed:
  n <- cond$n[pos]
  obs <- cond$obs[pos]
  lambda_noninvariance =  cond$lambda_noninvariance[pos] |> as.character()
  theta_noninvariance =  cond$theta_noninvariance[pos] |> as.character()
  nu_noninvariance =  cond$nu_noninvariance[pos] |> as.character()
  seed_cond <- cond$seed[pos]
  set.seed(seed_cond)
  
  #### set data generation parameters ####
  ## regression parameters:
  # group 1
  phi11_g1_pop <- .3
  phi22_g1_pop <- .3
  phi12_g1_pop <- .1
  phi21_g1_pop <- .1
  
  # group 2
  phi11_g2_pop <- .3
  phi22_g2_pop <- .3
  phi12_g2_pop <- .1
  phi21_g2_pop <- .1
  
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
  loadings_baseline <- .72
  
  # loadings in group 1 are all equal to baseline:
  lambda_g1 <- lavaan::lav_matrix_bdiag(list(rep(loadings_baseline, 4),
                                             rep(loadings_baseline, 4)))
  
  # loadings in group 2 are equal to group 1 unless modified by condition:
  lambda_g2 <- lambda_g1
  
  if(lambda_noninvariance == "uniform"){
    lambda_g2[3, 1] <- lambda_g2[7, 2] <- loadings_baseline*1.5
    lambda_g2[4, 1] <- lambda_g2[8, 2] <- loadings_baseline*1.5
  }
  if(lambda_noninvariance == "mixed"){
    lambda_g2[3, 1] <- lambda_g2[7, 2] <- loadings_baseline*1.5
    lambda_g2[4, 1] <- lambda_g2[8, 2] <- loadings_baseline*0.5
  }
  
  ## theta matrix (residual variances):
  resvar_baseline <- .48
  
  # residual variances in group 1 are all equal to baseline:
  theta_g1 <- matrix(0, nrow = 8, ncol = 8)
  diag(theta_g1) <- resvar_baseline
  
  # residual variances in group 2 are equal to group 1 unless modified by condition:
  theta_g2 <- theta_g1
  
  if(theta_noninvariance == "uniform"){
    theta_g2[3, 3] <- theta_g2[7, 7] <- resvar_baseline*0.5
    theta_g2[4, 4] <- theta_g2[8, 8] <- resvar_baseline*0.5
  }
  if(theta_noninvariance == "mixed"){
    theta_g2[3, 3] <- theta_g2[7, 7] <- resvar_baseline*0.5
    theta_g2[4, 4] <- theta_g2[8, 8] <- resvar_baseline*1.5
  }
  
  ## nu vector (intercepts)
  # intercepts are all 0 in group 1
  nu_g1 <- rep(0, 8)
  
  # intercepts in group 2 are equal to group 1 unless modified by condition:
  nu_g2 <- nu_g1
  
  if(nu_noninvariance == "uniform"){
    nu_g2[3] <- nu_g2[7] <- 1
    nu_g2[4] <- nu_g2[8] <- 1
  }
  if(nu_noninvariance == "mixed"){
    nu_g2[3] <- nu_g2[7] <- 1
    nu_g2[4] <- nu_g2[8] <- -1
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
      nu_i <- nu_g1
    }
    if(g == "group2"){
      phimat_i <- phimat_g2
      zetamat_i <- zetamat_g2
      lambda_i <- lambda_g2
      theta_i <- theta_g2
      nu_i <- nu_g2
    }
    
    # generate person-specific latent means:
    mu_i <- rnorm(2, mean = 5, sd = 1)
    
    # generate factor scores:
    eta_i <- sim_VAR(factors = 2, obs = obs,
                     phi = phimat_i, zeta = zetamat_i,
                     mu = mu_i,
                     burn_in = 20)
    
    # generate errors:
    epsilon_i <- mvrnorm(obs, mu = rep(0, 8),
                         Sigma = theta_i, empirical=T)
    
    nu_matrix <- rep(nu_i, each = obs) |> matrix(nrow = obs)
    
    # transform factor scores into observed scores:
    data_i <- as.matrix(eta_i[, c("eta1", "eta2")]) %*% t(lambda_i) + nu_matrix + epsilon_i |>
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
    "f1 =~ 0.72*v1 + v2 + v3 + v4
      v1 ~ 0*1
      f1 ~ NA*1",
    "f2 =~ 0.72*v5 + v6 + v7 + v8
      v5 ~ 0*1
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
  if(lambda_noninvariance == "no"){
    partial_noninvariances <- list(NULL, NULL)
  } else {
    partial_noninvariances <- list(c("f1 =~ v3", "f1 =~ v4"),
                                   c("f2 =~ v7", "f2 =~ v8"))
  }
  
  if(theta_noninvariance == "no"){
    partial_noninvariances <- partial_noninvariances
  } else {
    partial_noninvariances[[1]] <- c(partial_noninvariances[[1]], 
                                     "v3 ~~ v3", 
                                     "v4 ~~ v4")
    partial_noninvariances[[2]] <- c(partial_noninvariances[[2]], 
                                     "v7 ~~ v7", 
                                     "v8 ~~ v8")
  }
  
  if(nu_noninvariance == "no"){
    partial_noninvariances <- partial_noninvariances
  } else {
    partial_noninvariances[[1]] <- c(partial_noninvariances[[1]], 
                                     "v3 ~ 1", 
                                     "v4 ~ 1")
    partial_noninvariances[[2]] <- c(partial_noninvariances[[2]], 
                                     "v7 ~ 1", 
                                     "v8 ~ 1")
  }
  
  invariances <- c("loadings", "intercepts", "residuals")
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
    
    nu_ests <- output_step2_multi$result$result$MMparameters$nu_group
    nu_ests <- do.call(rbind, lapply(nu_ests, function(x){c(x[3:4, 1], x[7:8, 2])}))
    nu_true <- matrix(c(nu_g1[c(3, 4, 7, 8)], nu_g2[c(3, 4, 7, 8)]),
                      nrow = 2, byrow = TRUE)
    bias_nu <- sum(nu_true - nu_ests)/(4*2)  |> as.numeric()
    RMSE_nu <- sqrt(sum((nu_true - nu_ests)^2)/(4*2)) |> as.numeric()
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
    bias_nu <- NA
    RMSE_nu <- NA
  }
  
  duration <- difftime(Sys.time(), start, units = "secs") |> as.numeric()
  
  output <- c("iteration" = iteration, "replication" = replication,
              "n" = n, "obs" = obs,
              "lambda_noninvariance" = lambda_noninvariance,
              "theta_noninvariance" = theta_noninvariance,
              "nu_noninvariance" = nu_noninvariance,
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
              "bias_lambda" = bias_lambda, "bias_theta" = bias_theta, "bias_nu" = bias_nu,
              "RMSE_lambda" = RMSE_lambda, "RMSE_theta" = RMSE_theta, "RMSE_nu" = RMSE_nu,
              "step1_single_warning" = step1_single_warning, "step2_single_warning" = step2_single_warning, "step3_single_warning" = step3_single_warning,
              "step1_single_error" = step1_single_error, "step2_single_error" = step2_single_error, "step3_single_error" = step3_single_error,
              "step1_multi_warning" = step1_multi_warning, "step2_multi_warning" = step2_multi_warning, "step3_multi_warning" = step3_multi_warning,
              "step1_multi_error" = step1_multi_error, "step2_multi_error" = step2_multi_error, "step3_multi_error" = step3_multi_error,
              "seed" = seed_cond, "pos" = pos,
              "step1_single_warning_text" = step1_single_warning_text, "step2_single_warning_text" = step2_single_warning_text, "step3_single_warning_text" = step3_single_warning_text,
              "step1_single_error_text" = step1_single_error_text, "step2_single_error_text" = step2_single_error_text, "step3_single_error_text" = step3_single_error_text,
              "step1_multi_warning_text" = step1_multi_warning_text, "step2_multi_warning_text" = step2_multi_warning_text, "step3_multi_warning_text" = step3_multi_warning_text,
              "step1_multi_error_text" = step1_multi_error_text, "step2_multi_error_text" = step2_multi_error_text, "step3_multi_error_text" = step3_multi_error_text)
  
  for(i in 97:108){
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