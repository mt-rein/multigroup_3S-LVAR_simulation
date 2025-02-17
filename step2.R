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
  # Are there measurement blocks?
  if(is.list(measurementmodel)){
    for(m in 1:M){
      # if there are measurement blocks, compute factor scores in each block
      temp <- lavPredict(fit_step1[[m]], assemble = TRUE, append.data = TRUE) |> 
        as.data.frame()
      # and append to original data
      data <- dplyr::full_join(data, temp, by = c(group, indicators[indicators %in% colnames(temp)]))
    }
  } else {
    # if there are no measurement blocks, compute factor scores for all latent 
    # variables simultaneously and append to data
    temp <- lavPredict(fit_step1, assemble = TRUE, append.data = TRUE) |> 
      as.data.frame()
    data <- dplyr::full_join(data, temp, by = c(group, indicators))
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
    if(is.list(measurementmodel)){
      # if yes, create lists to store MM parameter values per block
      # (lists with M elements)
      psi_block <- lambda_block <- theta_block <- nu_block <- vector(mode = "list", length = M)
      for(m in 1:M){
        EST_block         <- lavaan::lavInspect(fit_step1[[m]], "est")
        psi_block[[m]] <- EST_block[["psi"]]
        lambda_block[[m]] <- EST_block[["lambda"]]
        theta_block[[m]]  <- EST_block[["theta"]]
        nu_block[[m]]  <- EST_block[["nu"]]
      }
      
      # combine the matrices of the different measurement blocks into a single matrix
      psi_group <- lavaan::lav_matrix_bdiag(psi_block)
      lambda_group <- lavaan::lav_matrix_bdiag(lambda_block)
      theta_group  <- lavaan::lav_matrix_bdiag(theta_block)
      nu_group  <- lavaan::lav_matrix_bdiag(nu_block)
      
      # name the matrices' rows and columns
      rownames(psi_group) <- colnames(psi_group) <- factors
      rownames(lambda_group) <- indicators
      colnames(lambda_group) <- factors
      rownames(theta_group) <- colnames(theta_group) <- indicators
      rownames(nu_group) <- indicators
      colnames(nu_group) <- factors
      
    } else {
      # if there are no measurement blocks, simply extract the MM parameters
      psi_group <- lambda_group <- theta_group <-  vector(mode = "list", length = G)
      psi_group <- lavaan::lavInspect(fit_step1, "est")$psi
      lambda_group <- lavaan::lavInspect(fit_step1, "est")$lambda
      theta_group <- lavaan::lavInspect(fit_step1, "est")$theta
      nu_group <- lavaan::lavInspect(fit_step1, "est")$nu
    }
    
    # compute lambda_star and theta_star
    sigma_g <- lambda_group %*% psi_group %*% t(lambda_group) + theta_group
    lambda_star[1, ] <- diag(psi_group %*% t(lambda_group) %*% solve(sigma_g) %*% lambda_group)
    theta_star[1, ] <- lambda_star[1, ]*(1-lambda_star[1, ])*diag(psi_group)
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
    
    # Are there measurement blocks?
    if(is.list(measurementmodel)){
      # if yes, create lists to store MM parameter values per block
      # (lists with M elements, where each element is a list with G elements)
      psi_block <- lambda_block <- theta_block <- nu_block <- vector(mode = "list", length = M)
      for(m in 1:M){
        EST_block         <- lavaan::lavInspect(fit_step1[[m]], "est")
        psi_block[[m]] <- lapply(X = EST_block, "[[", "psi")
        lambda_block[[m]] <- lapply(X = EST_block, "[[", "lambda")
        theta_block[[m]]  <- lapply(X = EST_block, "[[", "theta")
        nu_block[[m]]  <- lapply(X = EST_block, "[[", "nu")
      }
      
      # create lists to store MM parameter values per group
      # (i.e., combine the block-wise MM parameters per group)
      # --> Lists with G elements, where each element is a list with M elements
      psi_group <- lambda_group <- theta_group <- nu_group <- vector(mode = "list", length = G)
      for(g in 1:G){
        # put MM matrices of each group in the respective list
        for(m in 1:M) {
          psi_group[[g]][[m]] <- psi_block[[m]][[g]]
          lambda_group[[g]][[m]] <- lambda_block[[m]][[g]]
          theta_group[[g]][[m]]  <- theta_block[[m]][[g]]
          nu_group[[g]][[m]]  <- nu_block[[m]][[g]]
        }
        # combine the matrices of the different measurement blocks into a single matrix per group
        psi_group[[g]] <- lavaan::lav_matrix_bdiag(psi_group[[g]])
        lambda_group[[g]] <- lavaan::lav_matrix_bdiag(lambda_group[[g]])
        theta_group[[g]]  <- lavaan::lav_matrix_bdiag(theta_group[[g]])
        nu_group[[g]]  <- lavaan::lav_matrix_bdiag(nu_group[[g]])
        
        # name the matrices' rows and columns
        rownames(psi_group[[g]]) <- colnames(psi_group[[g]]) <- factors
        rownames(lambda_group[[g]]) <- indicators
        colnames(lambda_group[[g]]) <- factors
        rownames(theta_group[[g]]) <- colnames(theta_group[[g]]) <- indicators
        rownames(nu_group[[g]]) <- indicators
        colnames(nu_group[[g]]) <- factors
      }
    } else {
      # if there are no measurement blocks, simply extract the MM parameters per group
      psi_group <- lambda_group <- theta_group <-  vector(mode = "list", length = G)
      for(g in 1:G){
        psi_group[[g]] <- lavaan::lavInspect(fit_step1, "est")[[g]]$psi
        lambda_group[[g]] <- lavaan::lavInspect(fit_step1, "est")[[g]]$lambda
        theta_group[[g]] <- lavaan::lavInspect(fit_step1, "est")[[g]]$theta
        nu_group[[g]] <- lavaan::lavInspect(fit_step1, "est")[[g]]$nu
      }
    }
    
    # compute lambda_star and theta_star per group
    for(g in 1:G) {
      sigma_g <- lambda_group[[g]] %*% psi_group[[g]] %*% t(lambda_group[[g]]) + theta_group[[g]]
      A_g <- psi_group[[g]] %*% t(lambda_group[[g]]) %*% solve(sigma_g)
      lambda_star[g, ] <- diag(A_g %*% lambda_group[[g]])
      theta_star[g, ] <- diag(A_g %*% theta_group[[g]] %*% t(A_g))
    }
    
  }
  
  # assemble output
  output <- list("data" = data,
                 "lambda_star" = lambda_star,
                 "theta_star" = theta_star,
                 "other" = list("factors" = factors,
                                "indicators" =  indicators,
                                "step1group" = group),
                 "MMparameters" = list("psi_group" = psi_group,
                                       "lambda_group" = lambda_group,
                                       "theta_group" = theta_group,
                                       "nu_group" = nu_group))
  return(output)
}