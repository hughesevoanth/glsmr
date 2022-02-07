#' Test if a non-linear model is a better fit to the data than a linear model
#'
#' This function performs a linear and non-linear GAM and tests the two to determine if the non-linear model is a better fit to the data.
#'
#' @keywords GAM linear non-linear
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param outcome a single string character of the column name for the outcome or dependent or response variable
#' @param exposure a single string character of the column name for the exposure or independent or explanatory variable
#' @param linear_covariates a vector of string(s) that are also column names used to define variables that will be set as parametric (linear) covariates.
#' @param smooth_covariates a vector of string(s) that are also column names used to define variables that will be set as non-linear (smooth, s()) covariates.
#' @param rnt_outcome binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
#' @param weights_variable a single string character of the column name for a weights variable
#' @param outlier_cutoff a single numeric value to define a cutoff value for how many iqr or sd units outlier values
#' @param outlier_method a single string character of "iqr" or "sd" to determine if outlier should be determined by means and sd or medians and iqr.
#' @param messages should a progress message be printed to screen - binary TRUE or FALSE
#' @return returns a vector of summary statistics
#' @importFrom stats na.omit anova formula lm fitted.values
#' @export
#' @examples
#' non_linear_test()
non_linear_test = function( wdata,
                  outcome = NA,
                  exposure = NA,
                  linear_covariates = NA,
                  smooth_covariates = NA,
                  rnt_outcome = FALSE,
                  weights_variable = NA,
                  outlier_method = "iqr",
                  outlier_cutoff = 5,
                  messages = FALSE){

  ############################################
  ### 0. look for any errors in parameters
  ############################################
  # if (!strata[1] %in% c("quartiles","deciles") & length(strata) < 3 & class(strata) != "numeric" ){
  #   stop("strata parameter can either be defined as (1) quartiles, (2) deciles, or (3) a numeric vector of at least length 3 defining strata boundries")
  #  }

  ## PART I: SETTING UP THE DATA
  ############################################
  ### PART I: 1. Define Model Data
  ############################################
  # names(outcome) = "outcome"
  if(messages == TRUE){ message("1. Defining model data frame.") }

  model_variables = na.omit( c(outcome, exposure, linear_covariates, smooth_covariates, weights_variable) )
  mod_data = wdata[, c(model_variables)]

  ############################################
  ## PART I: 2. define covariates
  ############################################
  covariates = na.omit(c(linear_covariates, smooth_covariates))
  if(length(covariates) == 0){covariates = NULL}

  ############################################
  ## PART I: 3. Identify outcome outliers
  ############################################
  if(messages == TRUE){ message("2. Identifying outliers") }
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
  outliers = id_outliers( y = mod_data[, exposure], outlier_method = outlier_method, outlier_cutoff = outlier_cutoff)

  # How many exposure outliers
  number_of_exposure_outliers = length(outliers); names(number_of_exposure_outliers) = "number_of_exposure_outliers"

  ## Turn outliers to NA
  if(length(outliers) > 0){
    mod_data[outliers, exposure] = NA
  }

  ############################################
  ## PART I: 5. normality of outcome
  ############################################
  if(messages == TRUE){ message("3. Estimate Shapiro-Wilk normality W-stat") }
  W_outcome = normW(mod_data[, outcome]); names(W_outcome) = "W_outcome"

  ############################################
  ## PART I: 6. normality of exposure
  ############################################
  W_exposure = normW(mod_data[, exposure]); names(W_exposure) = "W_exposure"

  ############################################
  ## PART I: 7. rank normalize the outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    if(messages == TRUE){ message("4. Performing rank normal transformation of outcome") }
    mod_data[, outcome] = rntransform( mod_data[, outcome] )
  } else {
    if(messages == TRUE){ message("4. No rank normal transformation of outcome performed") }
  }

  ## PART II: OBSERVATIONAL MODELING
  ############################################
  ### PART II: 1. null GAM model
  ###             with no exposure smooth
  ############################################
  if(messages == TRUE){ message("7. running full observational NULL GAM model") }
  gam_mod0 = gamfit( wdata = mod_data,
                     dependent = outcome,
                     independent = NA,
                     linear_covariates = c(linear_covariates, exposure),
                     smooth_covariates = smooth_covariates)


  ############################################
  ### PART II: 2. full GAM model
  ############################################
  if(messages == TRUE){ message("6. running full observational GAM model") }
  gam_mod = gamfit( wdata = mod_data,
                    dependent = outcome,
                    independent = exposure,
                    linear_covariates = linear_covariates,
                    smooth_covariates = smooth_covariates)


  ############################################
  ## PART II: 3. ANOVA test of GAM with no
  ##             exposure smooth and full GAM
  ############################################
  if(messages == TRUE){ message("8. testing non-linearity of observational data: GAM vs NULL GAM") }
  a = anova(gam_mod0$fit, gam_mod$fit, test = "F")
  obs_nonlinearity_test = a[2,3:6]; names(obs_nonlinearity_test) = paste0( "", c("df","deviance","F","P"))

  ############################################
  ## PART II: 4. normality of residuals
  ############################################
  if(messages == TRUE){ message("9. Estimate Shapiro-Wilk normality W-stat of GAM residuals") }
  W_residuals = normW( residuals(gam_mod$fit) ); names(W_residuals) = "W_GAM_residuals"


  ############################################
  ### Pull data together
  ############################################
  out = unlist( c( unlist( obs_nonlinearity_test) ,
          number_of_outcome_outliers, number_of_exposure_outliers,
          W_outcome, W_exposure, W_residuals) )

  if(messages == TRUE){ message("9. returning results to user") }
  return(out)

} ## end of function
