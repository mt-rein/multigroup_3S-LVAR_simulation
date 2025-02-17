step1 <- function(data, measurementmodel, group = NULL,
                  invariances = NULL,
                  partial_noninvariances = NULL,
                  ...){
  # data:
  #   a data frame with the indicator and ID variables
  # measurementmodel:
  #   a string describing the measurement model using the lavaan syntax
  # id:
  #   a character that indicates the id variable (the variable that indicates
  #   which observations belong to which person)
  
  # Are there measurement blocks in the MM?
  if(is.list(measurementmodel)){
    # If there are measurement blocks, how many do we have?
    M <- length(measurementmodel)
    
    # create list, with elements equal to number of measurement blocks
    MMoutput <- vector("list", length = M)
    for(m in 1:M){
      # estimate MM in each block
      MMoutput[[m]] <- lavaan::cfa(measurementmodel[[m]],
                                   data = data,
                                   estimator = "ML",
                                   group = group,
                                   group.equal = invariances,
                                   group.partial = partial_noninvariances[[m]],
                                   se = "none",
                                   ...)
    }
  } else {
    # if there are no measurement blocks, estimate full MM
    MMoutput <- lavaan::cfa(measurementmodel,
                            data = data,
                            estimator = "ML",
                            group = group,
                            group.equal = invariances,
                            group.partial = partial_noninvariances,
                            se = "none",
                            ...)
  }
  
  # assemble output
  output <- list("MMoutput" = MMoutput,
                 "data" = data,
                 "measurementmodel" = measurementmodel)
  return(output)
}