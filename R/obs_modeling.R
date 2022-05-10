#' linear and non-linear (GAM) observational modeling
#'
#' This performs traditional linear and GAM or non-linear modeling
#'
#' @keywords GAM lm observational modeling
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
#' @return returns a obs_modeling vector of summary statistics
#' @importFrom stats na.omit anova quantile sd formula lm fitted.values quantile
#' @importFrom lmtest lrtest
#' @export
#' @examples
#' obs_modeling()
obs_modeling = function( wdata,
                         outcome = NA,
                         exposure = NA,
                         linear_covariates = NA,
                         smooth_covariates = NA,
                         rnt_outcome = FALSE,
                         weights_variable = NA,
                         outlier_method = "iqr",
                         outlier_cutoff = 5,
                         messages = FALSE){

  ## PART I: SETTING UP THE DATA
  if(messages == TRUE){ message("Part I. Setting up the data.") }

  ############################################
  ### PART I: 1. Define Model Data
  ############################################
  if(messages == TRUE){ message("Part I.1. Defining model data frame.") }
  model_variables = na.omit( c(outcome, exposure, linear_covariates, smooth_covariates, weights_variable) )
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
  ## PART I: 5. normality of outcome
  ############################################
  if(messages == TRUE){ message("Part I.5. Estimate Shapiro-Wilk normality W-stat for outcome") }
  W_outcome = normW(mod_data[, outcome]); names(W_outcome) = "W_outcome"

  ############################################
  ## PART I: 6. normality of exposure
  ############################################
  if(messages == TRUE){ message("Part I.6. Estimate Shapiro-Wilk normality W-stat for exposure") }
  W_exposure = normW(mod_data[, exposure]); names(W_exposure) = "W_exposure"

  ############################################
  ## PART I: 7. rank normalize the outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    if(messages == TRUE){ message("Part I.7. Performing rank normal transformation of outcome") }
    mod_data[, outcome] = rntransform( mod_data[, outcome] )
  } else {
    if(messages == TRUE){ message("Part I.7. No rank normal transformation of outcome performed") }
  }

  ## PART II: OBSERVATIONAL MODELING
  if(messages == TRUE){ message("Part II. Observational modeling") }

  ############################################
  ### PART II: 1. Linear model
  ############################################
  if(messages == TRUE){ message("Part II.1. lm() linear modeling") }
  lm_mod = lmfit( wdata = mod_data,
                  dependent = outcome,
                  independent = exposure,
                  covariates = covariates)
  names(lm_mod$summary) = paste0("lm_",names(lm_mod$summary))

  ## LINEAR MODELING WITH THE gam() function
  gam_mod_lm = gamfit( wdata = mod_data,
                        dependent = outcome,
                        independent = NA,
                        linear_covariates = c(smooth_covariates, linear_covariates, exposure),
                        smooth_covariates = NA)
  ## exposure estimates
  beta = summary(gam_mod_lm$fit)[1][[1]][exposure]; names(beta) = "exposure_beta"
  se = summary(gam_mod_lm$fit)[2][[1]][exposure]; names(se) = "exposure_se"
  tval = summary(gam_mod_lm$fit)[3][[1]][exposure]; names(tval) = "exposure_tval"
  p = summary(gam_mod_lm$fit)[4][[1]][exposure]; names(p) = "exposure_P"
  aic = AIC(gam_mod_lm$fit); names(aic) = "aic"
  loglik = logLik(gam_mod_lm$fit); names(loglik) = "loglik"
  gam_mod_lm$summary = c(gam_mod_lm$summary, beta, se, tval, p, aic, loglik)
  names(gam_mod_lm$summary) = paste0("gam_lm_", names(gam_mod_lm$summary))

  ############################################
  ### PART II: 1. null GAM model
  ###             with no exposure smooth
  ############################################
  if(messages == TRUE){ message("Part II.1. running full observational NULL GAM model") }
  gam_mod0 = gamfit( wdata = mod_data,
                     dependent = outcome,
                     independent = NA,
                     linear_covariates = c(linear_covariates, exposure),
                     smooth_covariates = smooth_covariates)
  ## Add AIC and LogLik to summary data
  loglik = logLik(gam_mod0$fit); names(loglik) = "loglik"
  aic = AIC(gam_mod0$fit); names(aic) = "aic"
  gam_mod0$summary = c(gam_mod0$summary, loglik, aic)
  ## edit names
  names(gam_mod0$summary) = paste0("gam0_",names(gam_mod0$summary))

  ############################################
  ### PART II: 2. full GAM model
  ############################################
  if(messages == TRUE){ message("Part II.2. running full observational GAM model") }
  gam_mod = gamfit( wdata = mod_data,
                    dependent = outcome,
                    independent = exposure,
                    linear_covariates = linear_covariates,
                    smooth_covariates = smooth_covariates)
  ## Add AIC and LogLik to summary data
  loglik = logLik(gam_mod$fit); names(loglik) = "loglik"
  aic = AIC(gam_mod$fit); names(aic) = "aic"
  gam_mod$summary = c(gam_mod$summary, loglik, aic)
  ## edit names
  names(gam_mod$summary) = paste0("gam_",names(gam_mod$summary))


  ############################################
  ## PART II: 3. ANOVA test of GAM with no
  ##             exposure smooth and full GAM
  ############################################
  if(messages == TRUE){ message("Part II.3. testing non-linearity of observational data: GAM vs NULL GAM") }
  ## F Test
  a = anova(gam_mod_lm$fit, gam_mod0$fit, gam_mod$fit, test = "F")
  Ftest = unlist( c(a[2,3:6], a[3,3:6]) ); names(Ftest) = c(paste0( "Ftest_lmVgam0_", c("df","deviance","F","P")),
                                                            paste0( "Ftest_gam0Vgam_", c("df","deviance","F","P")))
  ## LRT
  lrt = lmtest::lrtest( gam_mod_lm$fit, gam_mod0$fit, gam_mod$fit )
  lrt = unlist( c(lrt[2, 3:5], lrt[3, 3:5]) )
  names(lrt) = c( paste0( "lrt_lmVgam0_", c("df","ChisqStat","P") ),
                  paste0( "lrt_gam0Vgam_", c("df","ChisqStat","P") ) )



  if(messages == TRUE){ message("Part IV.2. returning results to user") }
  out = unlist( c(  lm_mod$summary, gam_mod_lm$summary, gam_mod0$summary, gam_mod$summary, Ftest, lrt) )

  return(out)

} ## end of function
