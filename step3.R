step3 <- function(step2output, id, step3group = NULL){
  # step2output:
  #   the object that was generated using the step2() function
  # ...

  #### 1) Preparations ####
  ## extract objects from step 1 output:
  data <- step2output$data
  lambda_star <- step2output$lambda_star
  theta_star <- step2output$theta_star
  factors <- step2output$other$factors
  factors_ind <- paste0(factors, "_ind")                                        # names of factor score variables (single indicators)
  step1group <- step2output$other$step1group                                    # name of the grouping variable in step1
  unique_ids <- unique(data[, id])                                              # vector of unique ids
  N <- length(unique_ids)

  #### 2) data manipulation ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data |> 
    dplyr::rename_with(~ factors_ind, all_of(factors))
  
  #### 3) create OpenMx matrices ####
  xdim <- length(factors)*2 # number of latent constructs in model is number of factors times 2 due to the random intercepts
  udim <- 1 # exogenous covariates (ignored so far)
  ydim <- length(factors) # number of indicators (i.e., factor score variables)
  
  ## A matrices (= dynamics)
  amat <- mxMatrix(type = "Full", nrow = xdim, ncol = xdim,
                   name = "A",
                   free = c(TRUE, TRUE, FALSE, FALSE,
                            TRUE, TRUE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE),
                   values = c(.1, .1, 0, 0,
                              .1, .1, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1),
                   labels = c("phi11", "phi12", NA, NA,
                              "phi21", "phi22", NA, NA,
                              NA, NA, NA, NA,
                              NA, NA, NA, NA),
                   lbound= c(-.9, -.9, NA, NA,
                             -.9, -.9, NA, NA,
                             NA, NA, NA, NA,
                             NA, NA, NA, NA),
                   ubound= c(.9, .9, NA, NA,
                             .9, .9, NA, NA,
                             NA, NA, NA, NA,
                             NA, NA, NA, NA),
                   byrow = TRUE)
  
  # B matrix (= exogenous covariates on latent constructs)
  bmat <- mxMatrix('Zero', nrow = xdim, ncol = udim,
                   name='B')
  
  # D matrix (= exogenous covariates on observed variables)
  dmat <- mxMatrix('Zero', nrow = ydim, ncol = udim,
                   name='D')
  
  # Q matrix (= innovation (co)variances)
  qmat <- mxMatrix("Symm", xdim, xdim,
                   name = "Q",
                   free = c(TRUE,
                            TRUE, TRUE,
                            FALSE, FALSE, FALSE,
                            FALSE, FALSE, FALSE, FALSE),
                   values = c(1,
                              .3, 1,
                              0, 0, 0,
                              0, 0, 0, 0),
                   labels = c("zeta1",
                              "zeta12", "zeta2",
                              NA, NA, NA,
                              NA, NA, NA, NA),
                   lbound= c(1e-6,
                             NA, 1e-6,
                             NA, NA, NA,
                             NA, NA, NA, NA),
                   byrow = TRUE)
  
  
  # x0 and P0 (= initial values and (co)variances of the latent constructs)
  xmat <- mxMatrix('Full', nrow = xdim, ncol = 1,
                   name='x0',
                   free = FALSE,
                   values = c(0, 0, 0, 0))
  
  pmat <- mxMatrix('Symm', nrow = xdim, ncol = xdim,
                   name='P0',
                   free = c(TRUE,
                            TRUE, TRUE,
                            FALSE, FALSE, TRUE,
                            FALSE, FALSE, FALSE, TRUE),
                   values = c(1e1,
                              1e1, 1e1,
                              0, 0, 1e1,
                              0, 0, 0, 1e1),
                   labels = c("P0_f1",
                              "P0_f1f2", "P0_f2",
                              NA, NA, "P0_icp1",
                              NA, NA, "P0_icp12", "P0_icp2"),
                   byrow = TRUE)
  
  # u (= covariates)
  umat <- mxMatrix('Zero', nrow = udim, ncol = 1,
                   name = 'u')
  
  #### 4) create OpenMx models ####
  ## if there is no grouping variable in step3:
  if(purrr::is_empty(step3group)){
    # create a list of models (one for each individual) for each latent class:
    personmodelnames <- paste0("id_", unique_ids)
    names(personmodelnames) <- unique_ids
    
    personmodel_list <- vector(mode = "list", length = N)
    for(i in unique_ids){
      # which step1 group (if any) does the individual belong to?
      if(!purrr::is_empty(step1group)){
        s1g <- data[data$id == i, step1group] |> unique()
      } else {
        s1g <- 1
      }
      
      # fix the loadings to lambda_star (potentially depending on step1 group)
      # C matrix
      cmat <- mxMatrix('Full', nrow = ydim, ncol = xdim,
                       name='C',
                       free = FALSE,
                       values = c(lambda_star[s1g, 1], 0, lambda_star[s1g, 1], 0,
                                  0, lambda_star[s1g, 2], 0, lambda_star[s1g, 2]),
                       labels = NA,
                       byrow = TRUE,
                       dimnames = list(factors_ind, c(paste0(factors),
                                                      paste0("intercept_", factors))
                       )
      )
      
      # fix the residual variances to theta_star (potentially depending on step 1 group)
      # R matrix (= measurement noise)
      rmat <- mxMatrix('Diag', nrow = ydim, ncol = ydim,
                       name = 'R',
                       free = FALSE,
                       values = theta_star[s1g, ],
                       labels = NA
      )
      
      # create the model
      modelname <- personmodelnames[i]  |> as.character()
      personmodel_list[[i]] <- mxModel(name = modelname,
                                       amat, bmat, cmat, dmat,
                                       qmat, rmat, xmat, pmat, 
                                       umat,
                                       mxExpectationStateSpace('A', 'B', 'C', 'D', 
                                                               'Q', 'R', 'x0', 'P0', 
                                                               'u'),
                                       mxFitFunctionML(),
                                       mxData(data[data$id == i, c(factors_ind)], 
                                              'raw'))
    }
    names(personmodel_list) <- personmodelnames
    
    # combine the person-models to a multi-subject model
    fullmodel <- mxModel("fullmodel",
                         personmodel_list,
                         mxFitFunctionMultigroup(personmodelnames))
    
    # generate starting values:
    fullmodel <- generate_startval(fullmodel)
    
    # fit the model
    fullmodelr <- mxRun(fullmodel)
  }
  
  ## if there is a grouping variable in step 3
  if(!purrr::is_empty(step3group)){
    # create a list of models (one for each individual) for each latent class:
    personmodelnames <- paste0("id_", unique_ids)
    names(personmodelnames) <- unique_ids
    
    personmodel_list <- vector(mode = "list", length = N)
    # create the person-models:
    for(g in unique(data[, step3group])){
      group_ids <- unique(data$id[data[, step3group] == g])
      for(i in group_ids){
        if(!purrr::is_empty(step1group)){
          s1g <- data[data$id == i, step1group] |> unique()
        } else {
          s1g <- 1
        }
        
        # C matrix (= factor loadings, here fixed to lambda_star)
        cmat <- mxMatrix('Full', nrow = ydim, ncol = xdim,
                         name='C',
                         free = FALSE,
                         values = c(lambda_star[s1g, 1], 0, lambda_star[s1g, 1], 0,
                                    0, lambda_star[s1g, 2], 0, lambda_star[s1g, 2]),
                         labels = NA,
                         byrow = TRUE,
                         dimnames = list(factors_ind, c(paste0(factors),
                                                        paste0("intercept_", factors))
                         )
        )
        
        # R matrix (= measurement noise, here fixed to theta_star)
        rmat <- mxMatrix('Diag', nrow = ydim, ncol = ydim,
                         name = 'R',
                         free = FALSE,
                         values = theta_star[s1g, ],
                         labels = NA
        )
        
        modelname <- personmodelnames[i]  |> as.character()
        temp_model <- mxModel(name = modelname,
                              amat, bmat, cmat, dmat,
                              qmat, rmat, xmat, pmat,
                              umat,
                              mxExpectationStateSpace('A', 'B', 'C', 'D',
                                                      'Q', 'R', 'x0', 'P0',
                                                      'u'),
                              mxFitFunctionML(),
                              mxData(data[data[, id] == i, c(factors_ind)],
                                     'raw'))
        temp_model <- omxSetParameters(temp_model,
                                       labels = names(coef(temp_model)),
                                       newlabels = paste0(names(coef(temp_model)),
                                                          "_", g))
        personmodel_list[[i]] <- temp_model
      }
    }
    names(personmodel_list) <- personmodelnames
    
    # combine the person-models to a multi-subject model
    fullmodel <- mxModel("fullmodel", personmodel_list,
                         mxFitFunctionMultigroup(personmodelnames))
    
    # generate starting values:
    fullmodel <- generate_startval(fullmodel)
    
    # fit the model
    fullmodelr <- mxTryHard(fullmodel)
  }
  
  #### 7) build the output ####
  estimates <- coef(fullmodelr)
  
  output <- list("data" = data,
                 "estimates" = estimates,
                 "model" = fullmodelr)
  return(output)
}