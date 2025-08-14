#### This script defines the simulation function do_sim() ####

do_sim2 <- function(pos, cond, outputfile, verbose = FALSE){
  # pos = position in the condition grid
  # cond = the condition grid
  # outputfile = file name for the output CSV file
  # verbose = if TRUE, prints a message after the iteration is finished
  
  #### for testing:
  # pos = 1
  # cond = cond_sim2
  # heterogeneity = "large"
  # nature = "continuous"
  
  replication <- cond$replication[pos]
  iteration <- cond$iteration[pos]
  # get condition levels and set seed:
  ss_t <- cond$ss_t[pos]
  heterogeneity = cond$heterogeneity[pos] |> as.character()
  nature = cond$nature[pos] |> as.character()
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
  
  SM_grouping <- c(rep("group1", 32), rep("group2", 32))
  if(nature == "categorical"){
    MM_grouping <- c(rep("subgroupA", 32*0.25), rep("subgroupB", 32*0.75), 
                     rep("subgroupA", 32*0.75), rep("subgroupB", 32*0.25))
  }
  
  ## baseline measurement parameters:
  loadings_g1 <- c(.75, .75, .75, .75)
  resvar_g1 <- c(.54, .34, .54, .34)
  tau_g1 <- c(1, 1, 1.5, 1.5)
  
  loadings_g2 <- c(.75, .75, .45, .45)
  resvar_g2 <- c(.74, .54, .74, .54)
  tau_g2 <- c(1, 1, 1, 1)
  
  
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
  
  true_MM <- data.frame(lambda3 = numeric(),
                        lambda4 = numeric(),
                        lambda7 = numeric(),
                        lambda8 = numeric(),
                        tau3 = numeric(),
                        tau4 = numeric(),
                        tau7 = numeric(),
                        tau8 = numeric(),
                        theta1 = numeric(),
                        theta2 = numeric(),
                        theta3 = numeric(),
                        theta4 = numeric(),
                        theta5 = numeric(),
                        theta6 = numeric(),
                        theta7 = numeric(),
                        theta8 = numeric())
  
  ## simulate the data:
  for(i in 1:64){
    ## generate factor scores:
    # get data-generating values from the correct group
    if(SM_grouping[i] == "group1"){
      phimat <- phimat_g1
      zetamat <- zetamat_g1
      mu <- c(3.09, .98)
      loadings <- loadings_g1
      resvar <- resvar_g1
      tau <- tau_g1
    }
    if(SM_grouping[i] == "group2"){
      phimat <- phimat_g2
      zetamat <- zetamat_g2
      mu <- c(4.56, .32)
      loadings <- loadings_g2
      resvar <- resvar_g2
      tau <- tau_g2
    }
    # generate person-specific latent means:
    mu_i <- c(rnorm(1, mu[1], 1),
              rnorm(1, mu[2], 1))
    
    # generate time series data (factor scores):
    eta_i <- sim_VAR(factors = 2, obs = ss_t,
                     phi = phimat, zeta = zetamat,
                     mu = mu_i,
                     burn_in = 20)
    
    ## generate observed scores:
    # introduce heterogeneity in the MM parameters:
    # if heterogeneity is continuous:
    if(nature == "continuous"){
      if(heterogeneity == "small"){
        loadingsf1_i <- loadings + c(0, 0, runif(2, -.1, .1))
        loadingsf2_i <- loadings + c(0, 0, runif(2, -.1, .1))
        resvarf1_i <- resvar + runif(4, -.1, .1)
        resvarf2_i <- resvar + runif(4, -.1, .1)
        tauf1_i <- tau + c(0, 0, runif(2, -.2, .2))
        tauf2_i <- tau + c(0, 0, runif(2, -.2, .2))
      }
      if(heterogeneity == "large"){
        loadingsf1_i <- loadings + c(0, 0, runif(2, -.2, .2))
        loadingsf2_i <- loadings + c(0, 0, runif(2, -.2, .2))
        resvarf1_i <- resvar + runif(4, -.2, .2)
        resvarf2_i <- resvar + runif(4, -.2, .2)
        tauf1_i <- tau + c(0, 0, runif(2, -.4, .4))
        tauf2_i <- tau + c(0, 0, runif(2, -.4, .4))
      }
    }
    
    # if heterogeneity is categorical:
    if(nature == "categorical"){
      if(MM_grouping[i] == "subgroupA"){
        if(heterogeneity == "small"){
          loadingsf1_i <- loadings + c(0, 0, .1, .1)
          loadingsf2_i <- loadings + c(0, 0, .1, .1)
          resvarf1_i <- resvar + .1
          resvarf2_i <- resvar + .1
          tauf1_i <- tau + c(0, 0, .2, .2)
          tauf2_i <- tau + c(0, 0, .2, .2)
        }
        if(heterogeneity == "large"){
          loadingsf1_i <- loadings + c(0, 0, .2, .2)
          loadingsf2_i <- loadings + c(0, 0, .2, .2)
          resvarf1_i <- resvar + .2
          resvarf2_i <- resvar + .2
          tauf1_i <- tau + c(0, 0, .4, .4)
          tauf2_i <- tau + c(0, 0, .4, .4)
        }
      }
      
      if(MM_grouping[i] == "subgroupB"){
        if(heterogeneity == "small"){
          loadingsf1_i <- loadings - c(0, 0, .1, .1)
          loadingsf2_i <- loadings - c(0, 0, .1, .1)
          resvarf1_i <- resvar - .1
          resvarf2_i <- resvar - .1
          tauf1_i <- tau - c(0, 0, .2, .2)
          tauf2_i <- tau - c(0, 0, .2, .2)
        }
        if(heterogeneity == "large"){
          loadingsf1_i <- loadings - c(0, 0, .2, .2)
          loadingsf2_i <- loadings - c(0, 0, .2, .2)
          resvarf1_i <- resvar - .2
          resvarf2_i <- resvar - .2
          tauf1_i <- tau - c(0, 0, .4, .4)
          tauf2_i <- tau - c(0, 0, .4, .4)
        }
      }
    }
    
    lambda_i <- lavaan::lav_matrix_bdiag(list(loadingsf1_i, loadingsf2_i))
    theta_i <- matrix(0, nrow = 8, ncol = 8)
    diag(theta_i) <- c(resvarf1_i, resvarf2_i)
    tau_i <- c(tauf1_i, tauf2_i)
    # generate errors (item residuals):
    epsilon_i <- mvrnorm(nrow(eta_i), mu = rep(0, 8),
                         Sigma = theta_i, empirical = FALSE)
    # transform factor scores into group's observed scores:
    data_i <- t(tau_i + lambda_i %*% t(eta_i[, c("eta1", "eta2")])) + epsilon_i |>
      as.data.frame()
    colnames(data_i) <- paste0("v", 1:8)
    data_i$id <- i
    data_i$obs <- eta_i$obs
    data_i$group <- SM_grouping[i]
    
    # merge group's observed items with full data frame:
    data <- dplyr::full_join(data, data_i, by = join_by(id, obs, group, v1, v2, v3, v4, v5, v6, v7, v8))
    
    ## extract data-generating MM parameters:
    true_MM[i,] <- c(loadingsf1_i[3:4], loadingsf2_i[3:4],
                     tauf1_i[3:4], tauf2_i[3:4],
                     resvarf1_i, resvarf2_i)
    
  }
  
  start <- Sys.time()
  model_step1 <- list(
    "f1 =~ 0.75*v1 + v2 + v3 + v4
      v1 ~ 1*1
      f1 ~ NA*1",
    "f2 =~ 0.75*v5 + v6 + v7 + v8
      v5 ~ 1*1
      f2 ~ NA*1")
  # which parameters are constrained to equality?
  invariances <- c("loadings", "intercepts")
  
  partial_noninvariances <- list(c("f1 =~ v3", "f1 =~ v4",
                                   "v3 ~ 1", "v4 ~ 1"),
                                 c("f2 =~ v7", "f2 =~ v8",
                                   "v7 ~ 1", "v8 ~ 1"))
  
  #### Person-specific ####
  output_step1_personspec <- run_step1(data = data,
                                      measurementmodel = model_step1,
                                      group = "id",
                                      invariances = invariances,
                                      partial_noninvariances = partial_noninvariances)
  
  # extract error/warning messages (if applicable):
  step1_personspec_warning <- ifelse(is_empty(output_step1_personspec$warnings),
                                    FALSE, TRUE)
  step1_personspec_warning_text <- ifelse(is_empty(output_step1_personspec$warnings),
                                         "",
                                         paste(c(output_step1_personspec$warnings),
                                               collapse = "; ")
  )
  step1_personspec_error <- ifelse(is_empty(output_step1_personspec$result$error),
                                  FALSE, TRUE)
  step1_personspec_error_text <- ifelse(is_empty(output_step1_personspec$result$error),
                                       "",
                                       paste(c(output_step1_personspec$result$error),
                                             collapse = "; "))
  
  #### Step 2 ####
  if(!step1_personspec_error){                                                             # only proceed if there is no error in step 1
    output_step2_personspec <- run_step2(step1output = output_step1_personspec$result$result)
    # extract error/warning messages (if applicable):
    step2_personspec_warning <- ifelse(is_empty(output_step2_personspec$warnings),
                                      FALSE, TRUE)
    step2_personspec_warning_text <- ifelse(is_empty(output_step2_personspec$warnings),
                                           "",
                                           paste(c(output_step2_personspec$warnings),
                                                 collapse = "; ")
    )
    step2_personspec_error <- ifelse(is_empty(output_step2_personspec$result$error),
                                    FALSE, TRUE)
    step2_personspec_error_text <- ifelse(is_empty(output_step2_personspec$result$error),
                                         "",
                                         paste(c(output_step2_personspec$result$error),
                                               collapse = "; ")
    )
  } else {
    step2_personspec_warning <- FALSE
    step2_personspec_warning_text <- "step1 not successful"
    step2_personspec_error <- FALSE
    step2_personspec_error_text <- "step1 not successful"
  }
  
  # # check whether the solution is admissible:
  if(grepl("variances are negative", step1_personspec_warning_text, ignore.case = TRUE) & !step2_personspec_error) {
    # if there are negative variances, remove those persons with negative variances and re-reun step 1 and step 2
    which_negative_theta <- which(sapply(output_step2_personspec$result$result$MMparameters$theta_group, function(mat) any(diag(mat) < 0)))
    which_negative_psi <- which(sapply(output_step2_personspec$result$result$MMparameters$psi_group, function(mat) any(diag(mat) < 0)))
    which_negative <- unique(c(which_negative_theta, which_negative_psi))
    data_subset <- data |> dplyr::filter(!id %in% which_negative)
    
    # rerun step 1 and step 2:
    output_step1_personspec <- run_step1(data = data_subset,
                                         measurementmodel = model_step1,
                                         group = "id",
                                         invariances = invariances,
                                         partial_noninvariances = partial_noninvariances)
    
    # extract error/warning messages (if applicable):
    step1_personspec_warning <- ifelse(is_empty(output_step1_personspec$warnings),
                                       FALSE, TRUE)
    step1_personspec_warning_text <- ifelse(is_empty(output_step1_personspec$warnings),
                                            "",
                                            paste(c(output_step1_personspec$warnings),
                                                  collapse = "; ")
    )
    step1_personspec_error <- ifelse(is_empty(output_step1_personspec$result$error),
                                     FALSE, TRUE)
    step1_personspec_error_text <- ifelse(is_empty(output_step1_personspec$result$error),
                                          "",
                                          paste(c(output_step1_personspec$result$error),
                                                collapse = "; "))
    
    if(!step1_personspec_error){                                                             # only proceed if there is no error in step 1
      output_step2_personspec <- run_step2(step1output = output_step1_personspec$result$result)
      # extract error/warning messages (if applicable):
      step2_personspec_warning <- ifelse(is_empty(output_step2_personspec$warnings),
                                         FALSE, TRUE)
      step2_personspec_warning_text <- ifelse(is_empty(output_step2_personspec$warnings),
                                              "",
                                              paste(c(output_step2_personspec$warnings),
                                                    collapse = "; ")
      )
      step2_personspec_error <- ifelse(is_empty(output_step2_personspec$result$error),
                                       FALSE, TRUE)
      step2_personspec_error_text <- ifelse(is_empty(output_step2_personspec$result$error),
                                            "",
                                            paste(c(output_step2_personspec$result$error),
                                                  collapse = "; ")
      )
    } else {
      step2_personspec_warning <- FALSE
      step2_personspec_warning_text <- "step1 not successful"
      step2_personspec_error <- FALSE
      step2_personspec_error_text <- "step1 not successful"
    }
    
    n_negvars_personspec <- length(which_negative)
    
  } else {
    n_negvars_personspec <- 0
  }
  
  #### Step 3 ####
  if(!step1_personspec_error & !step2_personspec_error){                                              # only proceed if there is no error in step 1 as well as step 2
    output_step3_personspec <- run_step3(step2output = output_step2_personspec$result$result,
                                        id = "id", step3group = "group")
    
    # extract error/warning messages (if applicable):
    step3_personspec_warning <- ifelse(is_empty(output_step3_personspec$warnings),
                                      FALSE, TRUE)
    step3_personspec_warning_text <- ifelse(is_empty(output_step3_personspec$warnings),
                                           "",
                                           paste(c(output_step3_personspec$warnings),
                                                 collapse = "; ")
    )
    step3_personspec_error <- ifelse(is_empty(output_step3_personspec$result$error),
                                    FALSE, TRUE)
    step3_personspec_error_text <- ifelse(is_empty(output_step3_personspec$result$error),
                                         "",
                                         paste(c(output_step3_personspec$result$error),
                                               collapse = "; ")
    )
    
    # check if the model converged:
    if(output_step3_personspec$result$result$model@output$status$code != 0){
      step3_personspec_error <- TRUE
      step3_personspec_error_text <- "step3 model estimation failed"
    }
    
    
  } else {
    step3_personspec_warning <- FALSE
    step3_personspec_warning_text <- "step1 or step2 not successful"
    step3_personspec_error <- FALSE
    step3_personspec_error_text <- "step1 or step2 not successful"
  }
  
  if(!step1_personspec_error & !step2_personspec_error & !step3_personspec_error){
    ## SM parameters:
    estimates <- output_step3_personspec$result$result$estimates
    standarderrors <- setNames(summary(output_step3_personspec$result$result$model)$parameters$Std.Error,
                               summary(output_step3_personspec$result$result$model)$parameters$name)
    
    personspec_phi11_g1 <- estimates["phi11_group1"] |> as.numeric()
    personspec_phi22_g1 <- estimates["phi22_group1"] |> as.numeric()
    personspec_phi12_g1 <- estimates["phi12_group1"] |> as.numeric()
    personspec_phi21_g1 <- estimates["phi21_group1"] |> as.numeric()
    
    personspec_zeta1_g1 <- estimates["zeta1_group1"] |> as.numeric()
    personspec_zeta2_g1 <- estimates["zeta2_group1"] |> as.numeric()
    personspec_zeta12_g1 <- estimates["zeta12_group1"] |> as.numeric()
    
    personspec_phi11_g2 <- estimates["phi11_group2"] |> as.numeric()
    personspec_phi22_g2 <- estimates["phi22_group2"] |> as.numeric()
    personspec_phi12_g2 <- estimates["phi12_group2"] |> as.numeric()
    personspec_phi21_g2 <- estimates["phi21_group2"] |> as.numeric()
    
    personspec_zeta1_g2 <- estimates["zeta1_group2"] |> as.numeric()
    personspec_zeta2_g2 <- estimates["zeta2_group2"] |> as.numeric()
    personspec_zeta12_g2 <- estimates["zeta12_group2"] |> as.numeric()
    
    personspec_phi11_g1_se <- standarderrors["phi11_group1"] |> as.numeric()
    personspec_phi22_g1_se <- standarderrors["phi22_group1"] |> as.numeric()
    personspec_phi12_g1_se <- standarderrors["phi12_group1"] |> as.numeric()
    personspec_phi21_g1_se <- standarderrors["phi21_group1"] |> as.numeric()
    
    personspec_zeta1_g1_se <- standarderrors["zeta1_group1"] |> as.numeric()
    personspec_zeta2_g1_se <- standarderrors["zeta2_group1"] |> as.numeric()
    personspec_zeta12_g1_se <- standarderrors["zeta12_group1"] |> as.numeric()
    
    personspec_phi11_g2_se <- standarderrors["phi11_group2"] |> as.numeric()
    personspec_phi22_g2_se <- standarderrors["phi22_group2"] |> as.numeric()
    personspec_phi12_g2_se <- standarderrors["phi12_group2"] |> as.numeric()
    personspec_phi21_g2_se <- standarderrors["phi21_group2"] |> as.numeric()
    
    personspec_zeta1_g2_se <- standarderrors["zeta1_group2"] |> as.numeric()
    personspec_zeta2_g2_se <- standarderrors["zeta2_group2"] |> as.numeric()
    personspec_zeta12_g2_se <- standarderrors["zeta12_group2"] |> as.numeric()
    
    
    ## MM parameters:
    lambda_ests <- output_step2_personspec$result$result$MMparameters$lambda_group
    lambda_ests <- do.call(rbind, lapply(lambda_ests, function(x){c(x[3:4, 1], x[7:8, 2])}))
    lambda_true <- true_MM[, 1:4]
    bias_lambda3 <- mean(lambda_ests[, 1] - lambda_true[, 1])
    bias_lambda4 <- mean(lambda_ests[, 2] - lambda_true[, 2])
    bias_lambda7 <- mean(lambda_ests[, 3] - lambda_true[, 3])
    bias_lambda8 <- mean(lambda_ests[, 4] - lambda_true[, 4])
    RMSE_lambda3 <- sqrt(mean((lambda_ests[, 1] - lambda_true[, 1])^2))
    RMSE_lambda4 <- sqrt(mean((lambda_ests[, 2] - lambda_true[, 2])^2))
    RMSE_lambda7 <- sqrt(mean((lambda_ests[, 3] - lambda_true[, 3])^2))
    RMSE_lambda8 <- sqrt(mean((lambda_ests[, 4] - lambda_true[, 4])^2))
    
    tau_ests <- output_step2_personspec$result$result$MMparameters$tau_group
    tau_ests <- do.call(rbind, lapply(tau_ests, function(x){c(x[3:4, 1], x[7:8, 2])}))
    tau_true <- true_MM[5:8]
    bias_tau3 <- mean(tau_ests[, 1] - tau_true[, 1])
    bias_tau4 <- mean(tau_ests[, 2] - tau_true[, 2])
    bias_tau7 <- mean(tau_ests[, 3] - tau_true[, 3])
    bias_tau8 <- mean(tau_ests[, 4] - tau_true[, 4])
    RMSE_tau3 <- sqrt(mean((tau_ests[, 1] - tau_true[, 1])^2))
    RMSE_tau4 <- sqrt(mean((tau_ests[, 2] - tau_true[, 2])^2))
    RMSE_tau7 <- sqrt(mean((tau_ests[, 3] - tau_true[, 3])^2))
    RMSE_tau8 <- sqrt(mean((tau_ests[, 4] - tau_true[, 4])^2))

    theta_ests <- output_step2_personspec$result$result$MMparameters$theta_group
    theta_ests <- do.call(rbind, lapply(theta_ests, diag))
    theta_true <- true_MM[9:16]
    bias_theta1 <- mean(theta_ests[, 1] - theta_true[, 1])
    bias_theta2 <- mean(theta_ests[, 2] - theta_true[, 2])
    bias_theta3 <- mean(theta_ests[, 3] - theta_true[, 3])
    bias_theta4 <- mean(theta_ests[, 4] - theta_true[, 4])
    bias_theta5 <- mean(theta_ests[, 5] - theta_true[, 5])
    bias_theta6 <- mean(theta_ests[, 6] - theta_true[, 6])
    bias_theta7 <- mean(theta_ests[, 7] - theta_true[, 7])
    bias_theta8 <- mean(theta_ests[, 8] - theta_true[, 8])
    RMSE_theta1 <- sqrt(mean((theta_ests[, 1] - theta_true[, 1])^2))
    RMSE_theta2 <- sqrt(mean((theta_ests[, 2] - theta_true[, 2])^2))
    RMSE_theta3 <- sqrt(mean((theta_ests[, 3] - theta_true[, 3])^2))
    RMSE_theta4 <- sqrt(mean((theta_ests[, 4] - theta_true[, 4])^2))
    RMSE_theta5 <- sqrt(mean((theta_ests[, 5] - theta_true[, 5])^2))
    RMSE_theta6 <- sqrt(mean((theta_ests[, 6] - theta_true[, 6])^2))
    RMSE_theta7 <- sqrt(mean((theta_ests[, 7] - theta_true[, 7])^2))
    RMSE_theta8 <- sqrt(mean((theta_ests[, 8] - theta_true[, 8])^2))
    
    
  } else {
    personspec_phi11_g1 <- NA
    personspec_phi22_g1 <- NA
    personspec_phi12_g1 <- NA
    personspec_phi21_g1 <- NA
    
    personspec_zeta1_g1 <- NA
    personspec_zeta2_g1 <- NA
    personspec_zeta12_g1 <- NA
    
    personspec_phi11_g2 <- NA
    personspec_phi22_g2 <- NA
    personspec_phi12_g2 <- NA
    personspec_phi21_g2 <- NA
    
    personspec_zeta1_g2 <- NA
    personspec_zeta2_g2 <- NA
    personspec_zeta12_g2 <- NA
    
    personspec_phi11_g1_se <- NA
    personspec_phi22_g1_se <- NA
    personspec_phi12_g1_se <- NA
    personspec_phi21_g1_se <- NA
    
    personspec_zeta1_g1_se <- NA
    personspec_zeta2_g1_se <- NA
    personspec_zeta12_g1_se <- NA
    
    personspec_phi11_g2_se <- NA
    personspec_phi22_g2_se <- NA
    personspec_phi12_g2_se <- NA
    personspec_phi21_g2_se <- NA
    
    personspec_zeta1_g2_se <- NA
    personspec_zeta2_g2_se <- NA
    personspec_zeta12_g2_se <- NA
    
    bias_lambda3 <- NA
    bias_lambda4 <- NA
    bias_lambda7 <- NA
    bias_lambda8 <- NA
    RMSE_lambda3 <- NA
    RMSE_lambda4 <- NA
    RMSE_lambda7 <- NA
    RMSE_lambda8 <- NA
    
    bias_tau3 <- NA
    bias_tau4 <- NA
    bias_tau7 <- NA
    bias_tau8 <- NA
    RMSE_tau3 <- NA
    RMSE_tau4 <- NA
    RMSE_tau7 <- NA
    RMSE_tau8 <- NA
    
    bias_theta1 <- NA
    bias_theta2 <- NA
    bias_theta3 <- NA
    bias_theta4 <- NA
    bias_theta5 <- NA
    bias_theta6 <- NA
    bias_theta7 <- NA
    bias_theta8 <- NA
    RMSE_theta1 <- NA
    RMSE_theta2 <- NA
    RMSE_theta3 <- NA
    RMSE_theta4 <- NA
    RMSE_theta5 <- NA
    RMSE_theta6 <- NA
    RMSE_theta7 <- NA
    RMSE_theta8 <- NA
  }
  
  
  #### Group-specific ####
  output_step1_groupspec <- run_step1(data = data,
                                  measurementmodel = model_step1,
                                  group = "group",
                                  invariances = invariances,
                                  partial_noninvariances = partial_noninvariances)
  
  # extract error/warning messages (if applicable):
  step1_groupspec_warning <- ifelse(is_empty(output_step1_groupspec$warnings),
                                FALSE, TRUE)
  step1_groupspec_warning_text <- ifelse(is_empty(output_step1_groupspec$warnings),
                                     "",
                                     paste(c(output_step1_groupspec$warnings),
                                           collapse = "; ")
  )
  step1_groupspec_error <- ifelse(is_empty(output_step1_groupspec$result$error),
                              FALSE, TRUE)
  step1_groupspec_error_text <- ifelse(is_empty(output_step1_groupspec$result$error),
                                   "",
                                   paste(c(output_step1_groupspec$result$error),
                                         collapse = "; "))
  
  # check whether the solution is admissible:
  if(grepl("variances are negative", step1_groupspec_warning_text, ignore.case = TRUE)) {
    # Re-run step 1 with wide bounds in case there's a heywood case
    output_step1_groupspec <- run_step1(data = data,
                                    measurementmodel = model_step1,
                                    group = "group",
                                    invariances = invariances,
                                    partial_noninvariances = partial_noninvariances,
                                    bounds = "wide")
    
    # extract error/warning messages (if applicable):
    step1_groupspec_warning <- ifelse(is_empty(output_step1_groupspec$warnings),
                                  FALSE, TRUE)
    step1_groupspec_warning_text <- ifelse(is_empty(output_step1_groupspec$warnings),
                                       "",
                                       paste(c(output_step1_groupspec$warnings),
                                             collapse = "; ")
    )
    step1_groupspec_error <- ifelse(is_empty(output_step1_groupspec$result$error),
                                FALSE, TRUE)
    step1_groupspec_error_text <- ifelse(is_empty(output_step1_groupspec$result$error),
                                     "",
                                     paste(c(output_step1_groupspec$result$error),
                                           collapse = "; "))
    rerun_step1_groupspec <- TRUE
  } else {
    rerun_step1_groupspec <- FALSE
  }
  
  #### Step 2 ####
  if(!step1_groupspec_error){                                                             # only proceed if there is no error in step 1
    output_step2_groupspec <- run_step2(step1output = output_step1_groupspec$result$result)
    # extract error/warning messages (if applicable):
    step2_groupspec_warning <- ifelse(is_empty(output_step2_groupspec$warnings),
                                  FALSE, TRUE)
    step2_groupspec_warning_text <- ifelse(is_empty(output_step2_groupspec$warnings),
                                       "",
                                       paste(c(output_step2_groupspec$warnings),
                                             collapse = "; ")
    )
    step2_groupspec_error <- ifelse(is_empty(output_step2_groupspec$result$error),
                                FALSE, TRUE)
    step2_groupspec_error_text <- ifelse(is_empty(output_step2_groupspec$result$error),
                                     "",
                                     paste(c(output_step2_groupspec$result$error),
                                           collapse = "; ")
    )
  } else {
    step2_groupspec_warning <- FALSE
    step2_groupspec_warning_text <- "step1 not successful"
    step2_groupspec_error <- FALSE
    step2_groupspec_error_text <- "step1 not successful"
  }
  
  #### Step 3 ####
  if(!step1_groupspec_error & !step2_groupspec_error){                                              # only proceed if there is no error in step 1 as well as step 2
    output_step3_groupspec <- run_step3(step2output = output_step2_groupspec$result$result,
                                    id = "id", step3group = "group")
    
    # extract error/warning messages (if applicable):
    step3_groupspec_warning <- ifelse(is_empty(output_step3_groupspec$warnings),
                                  FALSE, TRUE)
    step3_groupspec_warning_text <- ifelse(is_empty(output_step3_groupspec$warnings),
                                       "",
                                       paste(c(output_step3_groupspec$warnings),
                                             collapse = "; ")
    )
    step3_groupspec_error <- ifelse(is_empty(output_step3_groupspec$result$error),
                                FALSE, TRUE)
    step3_groupspec_error_text <- ifelse(is_empty(output_step3_groupspec$result$error),
                                     "",
                                     paste(c(output_step3_groupspec$result$error),
                                           collapse = "; ")
    )
    
    # check if the model converged:
    if(output_step3_groupspec$result$result$model@output$status$code != 0){
      step3_groupspec_error <- TRUE
      step3_groupspec_error_text <- "step3 model estimation failed"
    }
    
    
  } else {
    step3_groupspec_warning <- FALSE
    step3_groupspec_warning_text <- "step1 or step2 not successful"
    step3_groupspec_error <- FALSE
    step3_groupspec_error_text <- "step1 or step2 not successful"
  }
  
  if(!step1_groupspec_error & !step2_groupspec_error & !step3_groupspec_error){
    estimates <- output_step3_groupspec$result$result$estimates
    standarderrors <- setNames(summary(output_step3_groupspec$result$result$model)$parameters$Std.Error,
                               summary(output_step3_groupspec$result$result$model)$parameters$name)
    
    groupspec_phi11_g1 <- estimates["phi11_group1"] |> as.numeric()
    groupspec_phi22_g1 <- estimates["phi22_group1"] |> as.numeric()
    groupspec_phi12_g1 <- estimates["phi12_group1"] |> as.numeric()
    groupspec_phi21_g1 <- estimates["phi21_group1"] |> as.numeric()
    
    groupspec_zeta1_g1 <- estimates["zeta1_group1"] |> as.numeric()
    groupspec_zeta2_g1 <- estimates["zeta2_group1"] |> as.numeric()
    groupspec_zeta12_g1 <- estimates["zeta12_group1"] |> as.numeric()
    
    groupspec_phi11_g2 <- estimates["phi11_group2"] |> as.numeric()
    groupspec_phi22_g2 <- estimates["phi22_group2"] |> as.numeric()
    groupspec_phi12_g2 <- estimates["phi12_group2"] |> as.numeric()
    groupspec_phi21_g2 <- estimates["phi21_group2"] |> as.numeric()
    
    groupspec_zeta1_g2 <- estimates["zeta1_group2"] |> as.numeric()
    groupspec_zeta2_g2 <- estimates["zeta2_group2"] |> as.numeric()
    groupspec_zeta12_g2 <- estimates["zeta12_group2"] |> as.numeric()
    
    groupspec_phi11_g1_se <- standarderrors["phi11_group1"] |> as.numeric()
    groupspec_phi22_g1_se <- standarderrors["phi22_group1"] |> as.numeric()
    groupspec_phi12_g1_se <- standarderrors["phi12_group1"] |> as.numeric()
    groupspec_phi21_g1_se <- standarderrors["phi21_group1"] |> as.numeric()
    
    groupspec_zeta1_g1_se <- standarderrors["zeta1_group1"] |> as.numeric()
    groupspec_zeta2_g1_se <- standarderrors["zeta2_group1"] |> as.numeric()
    groupspec_zeta12_g1_se <- standarderrors["zeta12_group1"] |> as.numeric()
    
    groupspec_phi11_g2_se <- standarderrors["phi11_group2"] |> as.numeric()
    groupspec_phi22_g2_se <- standarderrors["phi22_group2"] |> as.numeric()
    groupspec_phi12_g2_se <- standarderrors["phi12_group2"] |> as.numeric()
    groupspec_phi21_g2_se <- standarderrors["phi21_group2"] |> as.numeric()
    
    groupspec_zeta1_g2_se <- standarderrors["zeta1_group2"] |> as.numeric()
    groupspec_zeta2_g2_se <- standarderrors["zeta2_group2"] |> as.numeric()
    groupspec_zeta12_g2_se <- standarderrors["zeta12_group2"] |> as.numeric()
  } else {
    groupspec_phi11_g1 <- NA
    groupspec_phi22_g1 <- NA
    groupspec_phi12_g1 <- NA
    groupspec_phi21_g1 <- NA
    
    groupspec_zeta1_g1 <- NA
    groupspec_zeta2_g1 <- NA
    groupspec_zeta12_g1 <- NA
    
    groupspec_phi11_g2 <- NA
    groupspec_phi22_g2 <- NA
    groupspec_phi12_g2 <- NA
    groupspec_phi21_g2 <- NA
    
    groupspec_zeta1_g2 <- NA
    groupspec_zeta2_g2 <- NA
    groupspec_zeta12_g2 <- NA
    
    groupspec_phi11_g1_se <- NA
    groupspec_phi22_g1_se <- NA
    groupspec_phi12_g1_se <- NA
    groupspec_phi21_g1_se <- NA
    
    groupspec_zeta1_g1_se <- NA
    groupspec_zeta2_g1_se <- NA
    groupspec_zeta12_g1_se <- NA
    
    groupspec_phi11_g2_se <- NA
    groupspec_phi22_g2_se <- NA
    groupspec_phi12_g2_se <- NA
    groupspec_phi21_g2_se <- NA
    
    groupspec_zeta1_g2_se <- NA
    groupspec_zeta2_g2_se <- NA
    groupspec_zeta12_g2_se <- NA
  }
  
  duration <- difftime(Sys.time(), start, units = "secs") |> as.numeric()
  
  output <- c("iteration" = iteration, "replication" = replication,
              "ss_t" = ss_t, "heterogeneity" = heterogeneity, "nature" = nature,
              "duration" = duration,
              "phi11_g1_pop" = phi11_g1_pop, "phi12_g1_pop" = phi12_g1_pop, "phi21_g1_pop" = phi21_g1_pop, "phi22_g1_pop" = phi22_g1_pop,
              "phi11_g2_pop" = phi11_g2_pop, "phi12_g2_pop" = phi12_g2_pop, "phi21_g2_pop" = phi21_g2_pop, "phi22_g2_pop" = phi22_g2_pop,
              "zeta1_g1_pop" = zeta1_g1_pop, "zeta2_g1_pop" = zeta2_g1_pop, "zeta12_g1_pop" = zeta12_g1_pop,
              "zeta1_g2_pop" = zeta1_g2_pop, "zeta2_g2_pop" = zeta2_g2_pop, "zeta12_g2_pop" = zeta12_g2_pop,
              "personspec_phi11_g1" = personspec_phi11_g1, "personspec_phi12_g1" = personspec_phi12_g1, "personspec_phi21_g1" = personspec_phi21_g1, "personspec_phi22_g1" = personspec_phi22_g1,
              "personspec_phi11_g2" = personspec_phi11_g2, "personspec_phi12_g2" = personspec_phi12_g2, "personspec_phi21_g2" = personspec_phi21_g2, "personspec_phi22_g2" = personspec_phi22_g2,
              "personspec_zeta1_g1" = personspec_zeta1_g1, "personspec_zeta2_g1" = personspec_zeta2_g1, "personspec_zeta12_g1" = personspec_zeta12_g1,
              "personspec_zeta1_g2" = personspec_zeta1_g2, "personspec_zeta2_g2" = personspec_zeta2_g2, "personspec_zeta12_g2" = personspec_zeta12_g2,
              "personspec_phi11_g1_se" = personspec_phi11_g1_se, "personspec_phi12_g1_se" = personspec_phi12_g1_se, "personspec_phi21_g1_se" = personspec_phi21_g1_se, "personspec_phi22_g1_se" = personspec_phi22_g1_se,
              "personspec_phi11_g2_se" = personspec_phi11_g2_se, "personspec_phi12_g2_se" = personspec_phi12_g2_se, "personspec_phi21_g2_se" = personspec_phi21_g2_se, "personspec_phi22_g2_se" = personspec_phi22_g2_se,
              "personspec_zeta1_g1_se" = personspec_zeta1_g1_se, "personspec_zeta2_g1_se" = personspec_zeta2_g1_se, "personspec_zeta12_g1_se" = personspec_zeta12_g1_se,
              "personspec_zeta1_g2_se" = personspec_zeta1_g2_se, "personspec_zeta2_g2_se" = personspec_zeta2_g2_se, "personspec_zeta12_g2_se" = personspec_zeta12_g2_se,
              "groupspec_phi11_g1" = groupspec_phi11_g1, "groupspec_phi12_g1" = groupspec_phi12_g1, "groupspec_phi21_g1" = groupspec_phi21_g1, "groupspec_phi22_g1" = groupspec_phi22_g1,
              "groupspec_phi11_g2" = groupspec_phi11_g2, "groupspec_phi12_g2" = groupspec_phi12_g2, "groupspec_phi21_g2" = groupspec_phi21_g2, "groupspec_phi22_g2" = groupspec_phi22_g2,
              "groupspec_zeta1_g1" = groupspec_zeta1_g1, "groupspec_zeta2_g1" = groupspec_zeta2_g1, "groupspec_zeta12_g1" = groupspec_zeta12_g1,
              "groupspec_zeta1_g2" = groupspec_zeta1_g2, "groupspec_zeta2_g2" = groupspec_zeta2_g2, "groupspec_zeta12_g2" = groupspec_zeta12_g2,
              "groupspec_phi11_g1_se" = groupspec_phi11_g1_se, "groupspec_phi12_g1_se" = groupspec_phi12_g1_se, "groupspec_phi21_g1_se" = groupspec_phi21_g1_se, "groupspec_phi22_g1_se" = groupspec_phi22_g1_se,
              "groupspec_phi11_g2_se" = groupspec_phi11_g2_se, "groupspec_phi12_g2_se" = groupspec_phi12_g2_se, "groupspec_phi21_g2_se" = groupspec_phi21_g2_se, "groupspec_phi22_g2_se" = groupspec_phi22_g2_se,
              "groupspec_zeta1_g1_se" = groupspec_zeta1_g1_se, "groupspec_zeta2_g1_se" = groupspec_zeta2_g1_se, "groupspec_zeta12_g1_se" = groupspec_zeta12_g1_se,
              "groupspec_zeta1_g2_se" = groupspec_zeta1_g2_se, "groupspec_zeta2_g2_se" = groupspec_zeta2_g2_se, "groupspec_zeta12_g2_se" = groupspec_zeta12_g2_se,
              "bias_lambda3" = bias_lambda3, "bias_lambda4" = bias_lambda4, "bias_lambda7" = bias_lambda7, "bias_lambda8" = bias_lambda8,
              "bias_tau3" = bias_tau3, "bias_tau4" = bias_tau4, "bias_tau7" = bias_tau7, "bias_tau8" = bias_tau8,
              "bias_theta1" = bias_theta1, "bias_theta2" = bias_theta2, "bias_theta3" = bias_theta3, "bias_theta4" = bias_theta4,
              "bias_theta5" = bias_theta5, "bias_theta6" = bias_theta6, "bias_theta7" = bias_theta7, "bias_theta8" = bias_theta8,
              "RMSE_lambda3" = RMSE_lambda3, "RMSE_lambda4" = RMSE_lambda4, "RMSE_lambda7" = RMSE_lambda7, "RMSE_lambda8" = RMSE_lambda8,
              "RMSE_tau3" = RMSE_tau3, "RMSE_tau4" = RMSE_tau4, "RMSE_tau7" = RMSE_tau7, "RMSE_tau8" = RMSE_tau8,
              "RMSE_theta1" = RMSE_theta1, "RMSE_theta2" = RMSE_theta2, "RMSE_theta3" = RMSE_theta3, "RMSE_theta4" = RMSE_theta4,
              "RMSE_theta5" = RMSE_theta5, "RMSE_theta6" = RMSE_theta6, "RMSE_theta7" = RMSE_theta7, "RMSE_theta8" = RMSE_theta8,
              "step1_personspec_warning" = step1_personspec_warning, "step2_personspec_warning" = step2_personspec_warning, "step3_personspec_warning" = step3_personspec_warning,
              "step1_personspec_error" = step1_personspec_error, "step2_personspec_error" = step2_personspec_error, "step3_personspec_error" = step3_personspec_error,
              "n_negvars_personspec" = n_negvars_personspec,
              "step1_groupspec_warning" = step1_groupspec_warning, "step2_groupspec_warning" = step2_groupspec_warning, "step3_groupspec_warning" = step3_groupspec_warning,
              "step1_groupspec_error" = step1_groupspec_error, "step2_groupspec_error" = step2_groupspec_error, "step3_groupspec_error" = step3_groupspec_error,
              "rerun_step1_groupspec" = rerun_step1_groupspec,
              "seed" = seed_cond, "pos" = pos,
              "step1_personspec_warning_text" = step1_personspec_warning_text, "step2_personspec_warning_text" = step2_personspec_warning_text, "step3_personspec_warning_text" = step3_personspec_warning_text,
              "step1_personspec_error_text" = step1_personspec_error_text, "step2_personspec_error_text" = step2_personspec_error_text, "step3_personspec_error_text" = step3_personspec_error_text,
              "step1_groupspec_warning_text" = step1_groupspec_warning_text, "step2_groupspec_warning_text" = step2_groupspec_warning_text, "step3_groupspec_warning_text" = step3_groupspec_warning_text,
              "step1_groupspec_error_text" = step1_groupspec_error_text, "step2_groupspec_error_text" = step2_groupspec_error_text, "step3_groupspec_error_text" = step3_groupspec_error_text)
  
  for(i in 125:136){
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