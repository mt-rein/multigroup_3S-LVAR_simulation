step2 <- function(step1output){
  # step1output:
  #   the object that was generated using the step1() function
  
  ## Preparations
  fit_step1 <- step1output$MMoutput
  data <- step1output$data
  measurementmodel <- step1output$measurementmodel
  indicators <- lavaan::lavNames(lavaan::lavaanify(measurementmodel), "ov")
  factors <- lavaan::lavNames(lavaan::lavaanify(measurementmodel), "lv")
  if(is.list(fit_step1)){
    group <- lavInspect(fit_step1[[1]], "group")
    M <- length(measurementmodel)
  } else {
    group <- lavInspect(fit_step1, "group")
  }
  
  if(!purrr::is_empty(group) & !is.character(data[, group])){
    #warning("Note: The step1 grouping variable has been transformed to a character.")
    data[, group] <- as.character(data[, group])
  }
  
  
  #### compute factor scores ####
  for(m in 1:M){
    # if there are measurement blocks, compute factor scores in each block
    temp <- lavPredict(fit_step1[[m]], assemble = TRUE, append.data = TRUE) |> 
      as.data.frame()
    # and append to original data
    data <- dplyr::full_join(data, temp, by = c(group, indicators[indicators %in% colnames(temp)]))
  }
  
  #### compute lambda_star and theta_star (no grouping variable) ####
  if(purrr::is_empty(group)){
    G <- 1
    
    # create empty matrices to store values.
    # One row per group (in this case, 1), one column per factor
    lambda_star <- theta_star <- matrix(NA, nrow = G, ncol = length(factors))
    # column names = names of the factors
    colnames(lambda_star) <- colnames(theta_star) <- factors
    
    # Are there measurement blocks?
    # create lists to store MM parameter values per block
    # (lists with M elements)
    psi_block <- alpha_block <- lambda_block <- theta_block <- tau_block <- vector(mode = "list", length = M)
    for(m in 1:M){
      EST_block         <- lavaan::lavInspect(fit_step1[[m]], "est")
      psi_block[[m]] <- EST_block[["psi"]]
      alpha_block[[m]] <- EST_block[["alpha"]]
      lambda_block[[m]] <- EST_block[["lambda"]]
      theta_block[[m]]  <- EST_block[["theta"]]
      tau_block[[m]]  <- EST_block[["nu"]]
    }
    
    # combine the matrices of the different measurement blocks into a single matrix
    psi_group <- lavaan::lav_matrix_bdiag(psi_block)
    alpha_group <- diag(length(factors))
    for(f in 1:length(factors)){
      alpha_group[f, f] <- alpha_block[[f]]
    }
    lambda_group <- lavaan::lav_matrix_bdiag(lambda_block)
    theta_group  <- lavaan::lav_matrix_bdiag(theta_block)
    tau_group  <- lavaan::lav_matrix_bdiag(tau_block)
    
    # name the matrices' rows and columns
    rownames(psi_group) <- colnames(psi_group) <- factors
    rownames(alpha_group) <- colnames(alpha_group) <- factors
    rownames(lambda_group) <- indicators
    colnames(lambda_group) <- factors
    rownames(theta_group) <- colnames(theta_group) <- indicators
    rownames(tau_group) <- indicators
    colnames(tau_group) <- factors
    
    # compute lambda_star and theta_star
    sigma_g <- lambda_group %*% psi_group %*% t(lambda_group) + theta_group
    A_g <- psi_group %*% t(lambda_group) %*% solve(sigma_g)
    lambda_star[1, ] <- diag(A_g %*% lambda_group)
    theta_star[1, ] <- diag(A_g %*% theta_group %*% t(A_g))
  }
  
  #### compute lambda_star and theta_star (with grouping variable) ####
  if(!purrr::is_empty(group)){
    # how many groups?
    G <- length(unique(data[, group]))
    
    # create empty matrices to store values.
    # One row per group, one column per factor
    lambda_star <- theta_star <- matrix(NA, nrow = G, ncol = length(factors))
    # add row/column names to the matrices
    # row-names = group labels
    # column names = names of the factors
    if(is.list(fit_step1)){
      rownames(lambda_star) <- rownames(theta_star) <- lavInspect(fit_step1[[1]], "group.label")
    } else {
      rownames(lambda_star) <- rownames(theta_star) <- lavInspect(fit_step1, "group.label")
    }
    colnames(lambda_star) <- colnames(theta_star) <- factors
    
    # create lists to store MM parameter values per block
    # (lists with M elements, where each element is a list with G elements)
    psi_block <- alpha_block <- lambda_block <- theta_block <- tau_block <- vector(mode = "list", length = M)
    for(m in 1:M){
      EST_block         <- lavaan::lavInspect(fit_step1[[m]], "est")
      psi_block[[m]] <- lapply(X = EST_block, "[[", "psi")
      alpha_block[[m]] <- lapply(X = EST_block, "[[", "alpha")
      lambda_block[[m]] <- lapply(X = EST_block, "[[", "lambda")
      theta_block[[m]]  <- lapply(X = EST_block, "[[", "theta")
      tau_block[[m]]  <- lapply(X = EST_block, "[[", "nu")
    }
    
    # create lists to store MM parameter values per group
    # (i.e., combine the block-wise MM parameters per group)
    # --> Lists with G elements, where each element is a list with M elements
    psi_group <- alpha_group <- lambda_group <- theta_group <- tau_group <- vector(mode = "list", length = G)
    for(g in 1:G){
      # put MM matrices of each group in the respective list
      for(m in 1:M) {
        psi_group[[g]][[m]] <- psi_block[[m]][[g]]
        alpha_group[[g]][[m]] <- alpha_block[[m]][[g]]
        lambda_group[[g]][[m]] <- lambda_block[[m]][[g]]
        theta_group[[g]][[m]]  <- theta_block[[m]][[g]]
        tau_group[[g]][[m]]  <- tau_block[[m]][[g]]
      }
      # combine the matrices of the different measurement blocks into a single matrix per group
      psi_group[[g]] <- lavaan::lav_matrix_bdiag(psi_group[[g]])
      alpha_group[[g]] <- lavaan::lav_matrix_bdiag(alpha_group[[g]])
      lambda_group[[g]] <- lavaan::lav_matrix_bdiag(lambda_group[[g]])
      theta_group[[g]]  <- lavaan::lav_matrix_bdiag(theta_group[[g]])
      tau_group[[g]]  <- lavaan::lav_matrix_bdiag(tau_group[[g]])
      
      # name the matrices' rows and columns
      rownames(psi_group[[g]]) <- colnames(psi_group[[g]]) <- factors
      rownames(alpha_group[[g]]) <- colnames(alpha_group[[g]]) <- factors
      rownames(lambda_group[[g]]) <- indicators
      colnames(lambda_group[[g]]) <- factors
      rownames(theta_group[[g]]) <- colnames(theta_group[[g]]) <- indicators
      rownames(tau_group[[g]]) <- indicators
      colnames(tau_group[[g]]) <- factors
    }
    
    
    # compute lambda_star and theta_star per group
    for(g in 1:G) {
      sigma_g <- lambda_group[[g]] %*% psi_group[[g]] %*% t(lambda_group[[g]]) + theta_group[[g]]
      A_g <- psi_group[[g]] %*% t(lambda_group[[g]]) %*% solve(sigma_g)
      lambda_star[g, ] <- diag(A_g %*% lambda_group[[g]])
      theta_star[g, ] <- diag(A_g %*% theta_group[[g]] %*% t(A_g))
    }
    
  }
  
  MMparameters <- list("psi_group" = psi_group,
                       "alpha_group" = alpha_group,
                       "lambda_group" = lambda_group,
                       "theta_group" = theta_group,
                       "tau_group" = tau_group)
  
  # assemble output
  output <- list("data" = data,
                 "lambda_star" = lambda_star,
                 "theta_star" = theta_star,
                 "other" = list("factors" = factors,
                                "indicators" =  indicators,
                                "step1group" = group),
                 "MMparameters" = MMparameters)
  return(output)
}