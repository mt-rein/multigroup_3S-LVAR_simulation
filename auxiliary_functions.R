#### This script defines auxiliary functions for the simulation ####

#### sim_VAR() ####
# this function generates data for a single individual according to a vector autoregressive model
sim_VAR <- function(factors, obs, phi, zeta, mu, burn_in = 0){
  # factors = number of factors
  # obs = number of observations
  # phi = auto-regressive effect (a matrix in case of multiple constructs)
  # zeta = innovation variance (a matrix in case of multiple constructs)
  # mu = latent means (a vector in case of multiple constructs)
  # burn_in = length of burn in (i.e., data that are generated to remove influence of initial random draw)
  
  
  # create empty dataframe of length obs + burn_in
  data <- as.data.frame(matrix(NA, nrow = burn_in + obs, ncol = factors))
  names(data) <- paste0("eta", 1:factors)
  
  intercept <- (diag(factors) - phi) %*% mu
  
  totalvar <- solve(diag(factors*factors) - kronecker(phi,phi)) %*% c(zeta) |> matrix(ncol = factors)
  
  for(i in 1:nrow(data)){
    # simulate the first observation from the person's starting point (= intercept)
    if(i == 1){
      data[i,] <- MASS::mvrnorm(n = 1, mu = mu, Sigma = totalvar)
    }
    
    # then loop through all the rows, predict the current observation from the previous observations, then add random innovation
    if(i > 1){
      predicted <- intercept + phi %*% unlist(data[i-1,])
      data[i, ] <- MASS::mvrnorm(n = 1, mu = predicted, Sigma = zeta)
    }
  }
  
  # remove the first rows, depending on length of burn in
  if(burn_in > 0){
    data <- dplyr::slice(data, burn_in+1:n()) 
  }
  
  data$obs <- 1:nrow(data)
  
  return(data)
}

#### generate_startval() ####
# generate starting values
generate_startval <- function(model){
  values <- coef(model)
  values[grep("^phi", names(values))] <- runif(length(grep("^phi", names(values))), .05, .3)
  values[grep("zeta[0-9]_+", names(values))] <- runif(length(grep("zeta[0-9]_+", names(values))), .5, 1.5)
  values[grep("zeta12_+", names(values))] <- runif(length(grep("zeta12_+", names(values))), .1, .5)
  model <- omxSetParameters(model,
                            labels = names(values),
                            values = values)
  return(model)
}

#### safely/quietly functions ####
run_step1 <- quietly(safely(step1))
run_step2 <- quietly(safely(step2))
run_step3 <- quietly(safely(step3))