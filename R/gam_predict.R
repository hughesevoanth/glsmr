#' Generalized additive model fitting
#'
#' This function fits a GAM to the data
#'
#' @keywords GAM
#' @param wdata a data frame passed to function containing necessary data for analysis.
#' @param outcome a single string character of the column name for the dependent or response variable.
#' @param exposure a single string character of the column name for the independent or explanatory variable.
#' @param linear_covariates a vector of string(s) that are also column names used to define variables that will be set as parametric (linear) covariates.
#' @param smooth_covariates a vector of string(s) that are also column names used to define variables that will be set as non-linear (smooth, s()) covariates.
#' @param rnt_dependent binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
#' @param new_pred_data_frame a data matrix containing the model independent and covariates from which exposure values should be predicted. If left as NA, an equal distant range of dependent values of length 25, and the mean of each numeric and max of each factor will be used in the prediction matrix.
#' @param prediction_n how many values should be predicted across the range of the exposure? Default is 25.
#' @param return_prediction_only binary TRUE or FALSE, defaulted to TRUE, asking if only a vector of exposure predicted values should be returned or a data frame of predicted values and the variable data used in the prediction.
#' @param messages should verbose messages be printed to screen
#' @return returns a list object containing two objects, (1) 'fit' which holds the gam() model and (2) 'summary' which is a vector of summary statistics derived from the model.
#' @importFrom mgcv gam s
#' @importFrom stats na.omit formula residuals
#' @export
#' @examples
#' gam_predict()
gam_predict = function( wdata,
                  outcome = NA,
                  exposure = NA,
                  linear_covariates = NA,
                  smooth_covariates = NA,
                  rnt_outcome = FALSE,
                  return_prediction_only = TRUE,
                  prediction_n = 25,
                  outlier_cutoff = 5,
                  outlier_method = "iqr",
                  new_pred_data_frame = NA,
                  messages = FALSE){

  ## PART I: SETTING UP THE DATA
  if(messages == TRUE){ message("Part I. Setting up the data.") }

  ############################################
  ### PART I: 1. Define Model Data
  ############################################
  if(messages == TRUE){ message("Part I.1. Defining model data frame.") }
  model_variables = na.omit( c(outcome, exposure, linear_covariates, smooth_covariates) )
  mod_data =  wdata[, c(model_variables)]

  ############################################
  ## PART I: 2. define covariates
  ############################################
  if(messages == TRUE){ message("Part I.2. Defining covariates.") }
  covariates = na.omit(c(linear_covariates, smooth_covariates))
  if(length(covariates) == 0){covariates = NULL}

  ############################################
  ## PART I: 3. Identify outcome outliers
  ############################################
  if(messages == TRUE){ message("Part I.3. Identifying outcome outliers") }
  outliers = id_outliers( y = mod_data[, outcome], outlier_method = outlier_method, outlier_cutoff = outlier_cutoff)

  # How many outcome outliers
  number_of_outcome_outliers = length(outliers); names(number_of_outcome_outliers) = "number_of_outcome_outliers"

  ## Turn outliers to NA
  if(length(outliers) > 0){
    mod_data[outliers, outcome] = NA
  }

  ############################################
  ## PART I: 4. Identify exposure outliers
  ############################################
  if(messages == TRUE){ message("Part I.4. Identifying exposure outliers") }
  outliers = id_outliers( y = mod_data[, exposure], outlier_method = outlier_method, outlier_cutoff = outlier_cutoff)

  # How many exposure outliers
  number_of_exposure_outliers = length(outliers); names(number_of_exposure_outliers) = "number_of_exposure_outliers"

  ## Turn outliers to NA
  if(length(outliers) > 0){
    mod_data[outliers, exposure] = NA
  }

  ############################################
  ## PART I: 5. rank normalize the outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    if(messages == TRUE){ message("Part I.5. Performing rank normal transformation of outcome") }
    mod_data[, outcome] = rntransform( mod_data[, outcome] )
  } else {
    if(messages == TRUE){ message("Part I.5. No rank normal transformation of outcome performed") }
  }

  ## PART II: OBSERVATIONAL MODELING
  if(messages == TRUE){ message("Part II. Fitting the GAM") }
  ############################################
  ### PART II: 2. full GAM model
  ############################################
  if(messages == TRUE){ message("Part II.2. Running full observational GAM model") }
  gam_mod = gamfit( wdata = mod_data,
                    dependent = outcome,
                    independent = exposure,
                    linear_covariates = linear_covariates,
                    smooth_covariates = smooth_covariates)


  ############################################
  ## PART III: Outcome Prediction
  ############################################
  if(messages == TRUE){ message("Part III. Outcome Prediction") }

  ###############################
  ## IF the user passes their own
  ## data frame of values to predict from
  ###############################
  if( !is.na(new_pred_data_frame) ){
    predicted_tbl = mgcv::predict.gam( gam_mod$fit, new_data, se.fit = TRUE)
  }

  ###############################
  ## IF the user does NOT passe their
  ## own data frame of values to predict from
  ###############################
  if( is.na(new_pred_data_frame) ){
    ## tidymv prediction FUNCTION
    predicted_tbl = tidymv::get_gam_predictions( gam_mod$fit,
                                                 series = {{exposure}} ,
                                                 series_length = prediction_n )
    ## limit to one value of each linear covariate
    if(!is.na(linear_covariates)){
      for(lc in linear_covariates){
        w = which(colnames(predicted_tbl) == lc)
        vals = unique(predicted_tbl[,w])
        if(length(vals) == 2){ ## if length is 2 as it might be for sex, just take the first one.
          vals = vals[1]
        } else { ## take the middle value
          vals = vals[ round( mean( 1:length(vals) ) ) ]
        }
        ## Reduce
        predicted_tbl = predicted_tbl %>% filter( .data[[lc]] == vals )
      }
    }
  }

  ####################
  ## Return results
  ####################
  if(return_prediction_only == TRUE){
    return(predicted_tbl[, outcome])
  } else {
    return(predicted_tbl)
  }

}

   # ############################################
  # ### PART III: Prediction
  # ############################################
  # if(messages == TRUE){ message("Part III.1. Checking for a data frame of new data in variable new_pred_data_frame") }
  #
  # if( is.na(new_pred_data_frame) ){
  #   if(messages == TRUE){ message( paste0("\t.... a data frame not provided. That's OK.") ) }
  #   if(messages == TRUE){ message("Part III.2. Data to perform prediction will be auto generated.") }
  #
  #   ## From the covariates identify what is a factor and what is not
  #   temp = mod_data[, covariates]
  #   cl = sapply(temp, class)
  #
  #   ## Are any characters?
  #   if(messages == TRUE){ message( paste0("\t.... identifying the most common string in each factor present in data set.") ) }
  #   common_factor_values = c()
  #   f = sapply(cl, function(x){ x == "character" | x == "factor"})
  #   if(sum(f)>0){
  #     w = which(f == TRUE)
  #     common_factor_values = sapply(w, function(i){
  #       counts = table( temp[,i] )
  #       value = names(which(counts == max(counts) ))
  #       return(value)
  #     })
  #   }
  #
  #   ## Are any numeric?
  #   if(messages == TRUE){ message( paste0("\t.... identifying the most common value (the mean) of each numeric|integer present in data set.") ) }
  #   common_numeric_values = c()
  #   f = sapply(cl, function(x){ x == "numeric" | x == "integer"})
  #   if(sum(f)>0){
  #     w = which(f == TRUE)
  #     common_numeric_values = sapply(w, function(i){
  #       value = mean( temp[,i], na.rm = TRUE )
  #       return(value)
  #     })
  #   }
  #
  #   ## define new data to predict
  #   new_data = mod_data[, c(covariates, exposure)]
  #
  #
  #   ## what is the range of exposure to predict?
  #   if(messages == TRUE){ message( paste0("\t.... identifying the the range and sequence of values for the exposure.") ) }
  #   range = quantile( new_data[, exposure] , probs = c(0, 1), na.rm = TRUE)
  #   exposure_vals = seq(range[1], range[2], by = (range[2]-range[1])/(prediction_n-1) )
  #
  #   ## RE-define new_data
  #   new_data = new_data[ 1:prediction_n, ]
  #   new_data[, exposure] = exposure_vals
  #
  #   ## redefine numerics
  #   if(length(common_numeric_values) > 0 ){
  #     for(i in 1:length(common_numeric_values)){
  #       w = which(colnames(new_data) == names(common_numeric_values)[i])
  #       new_data[,w] = common_numeric_values[i]
  #     }
  #   }
  #
  #   ## redefine factors
  #   if(length(common_factor_values) > 0 ){
  #     for(i in 1:length(common_factor_values)){
  #       w = which(colnames(new_data) == names(common_factor_values)[i])
  #       new_data[,w] = common_factor_values[i]
  #     }
  #   }
  # } else {
  #   if(messages == TRUE){ message("Part III.2. Outcome prediction will be generated using data frame passed to function.") }
  #   new_data = new_pred_data_frame
  # }
  #
  # ## Prediction
  # if(messages == TRUE){ message("Part III.3. Performing prediction.") }
  # pred <- mgcv::predict.gam( gam_mod$fit, new_data, se.fit = TRUE)
  # out = data.frame(predicted_outcome = pred,  new_data)
  #
  # ## Return prediction to user
  # if(messages == TRUE){ message("Part IV. Returning predictions to user.") }
  #
  # if(return_prediction_only == TRUE){
  #   return( out$predicted_outcome )
  # } else {
  #   return(out)
  # }
  #
# }
