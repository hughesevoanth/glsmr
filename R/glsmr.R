#' GAM and linear-strata Mendelian randomization
#'
#' This function derives strata defined observational and TSLS (MR) estimates using linear modeling.
#' The full data set will be fit to a linear model and a generalized linear model (GAM).
#' Then subsequently the exposure data will be stratified and a linear model will be used to derive
#' effect estimates for each strata.
#'
#' @keywords GAM strata stratified MR
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param outcome a single string character of the column name for the outcome or dependent or response variable
#' @param exposure a single string character of the column name for the exposure or independent or explanatory variable
#' @param instrument a data frame passed to function containing necessary data for analysis
#' @param linear_covariates a vector of string(s) that are also column names used to define variables that will be set as parametric (linear) covariates.
#' @param smooth_covariates a vector of string(s) that are also column names used to define variables that will be set as non-linear (smooth, s()) covariates.
#' @param strata a string or vector to define how strata should be defined.
#' * `strata = quartiles`: automatically identifies quartiles or 4 equally sized strata
#' * `strata = decile`: automatically identifies deciles or 10 equally sized strata
#' * `strata = c(1,10,20,30)`: a user defined numeric vector to define boundaries for each strata.
#'    - The numeric vector example provided will define 4 strata. Lower bound values are inclusive, upper bounds are exclusive, to the exception of the last bounding value.
#' @param rnt_outcome binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
#' @param weights_variable a single string character of the column name for a weights variable
#' @param outlier_cutoff a single numeric value to define a cutoff value for how many iqr or sd units outlier values
#' @param outlier_method a single string character of "iqr" or "sd" to determine if outlier should be determined by means and sd or medians and iqr.
#' @param messages should a progress message be printed to screen - binary TRUE or FALSE
#' @return returns a glsmr object containing the complete linear and GAM models for the full data set, summary statistics for the data, a strata observational table, and a strata TSLS (MR) table.
#' @importFrom stats na.omit anova quantile sd formula lm fitted.values quantile
#' @export
#' @examples
#' glsmr()
glsmr = function( wdata,
                  outcome = NA,
                  exposure = NA,
                  instrument = NA,
                  linear_covariates = NA,
                  smooth_covariates = NA,
                  strata = 4,
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

  model_variables = na.omit( c(outcome, exposure, instrument, linear_covariates, smooth_covariates, weights_variable) )
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
  ## PART II: 4. exposure summary statistics
  ############################################
  temp = na.omit( mod_data[, c(outcome, exposure, covariates)] )
  exp_ss = c( N = nrow(temp) ,
              mean = mean(temp[, exposure], na.rm = TRUE),
              min = min(temp[, exposure], na.rm = TRUE),
              max = max(temp[, exposure], na.rm = TRUE),
              sd = sd(temp[, exposure], na.rm = TRUE) )
  rm(temp)

  ############################################
  ### PART II: 5. Linear model
  ############################################
  if(messages == TRUE){ message("5. running full observational linear model") }
  lm_mod = lmfit( wdata = mod_data,
                  dependent = outcome,
                  independent = exposure,
                  covariates = covariates)

  ############################################
  ### PART II: 6. Combine summary stats of full
  ###             for data observational data
  ############################################
  obs_ss = c(lm_mod$summary, exp_ss)

  ############################################
  ## PART II: 7. Stratify data by Exposure
  ############################################
  if(messages == TRUE){ message("9. stratifying the observational data by exposure") }
  strata_data = stratify_data( wdata = mod_data, stratify_on = exposure, strata = strata )

  ############################################
  ## PART II: 8. Summary Statistics for exposure
  ##             by strata
  ############################################
  if(messages == TRUE){ message("10. stratified summary statistic for exposure") }
  strata_exp_ss = stratify_sumstats(wdata = strata_data, exposure = exposure)

  ############################################
  ## PART II: 9. Run observational linear model on each
  ##             strata
  ############################################
  if(messages == TRUE){ message("11. running observational linear models on each strata") }
  strata_lm_mod = stratify_lmfit( wdata = strata_data,
                               outcome = outcome,
                               exposure = exposure,
                               covariates = covariates)


  ############################################
  ## PART II: 10. Combine strata sumstats and
  ##              linear model estimates
  ############################################
  strata_obs_ss = cbind(strata_lm_mod, strata_exp_ss)

  ############################################
  ## PART II: 11. Combine strata summary stats
  ##              & full model summary stats
  ############################################
  obs_ss = rbind( strata_obs_ss , fulldata = obs_ss )
  obs_ss = obs_ss[, -c(4:6,9) ]

  ## PART III: Two stage least square | MR modeling
  ############################################
  ### PART III: 1. Derive Instrument free exposure
  ############################################
  if(messages == TRUE){ message("9. derive instument free exposure") }
  mod_data = iv_free_exposure( wdata = mod_data,
                    exposure = exposure,
                    instrument = instrument,
                    covariates = covariates,
                    exposure_mean_normalize = TRUE)

  ####################################
  ## PART III: 2. Derive Instrument predicted exposure or d.hat
  ####################################
  if(messages == TRUE){ message("9. derive instument predicted exposure") }
  mod_data = iv_predicted_exposure( wdata = mod_data,
                               exposure = exposure,
                               instrument = instrument,
                               covariates = covariates)

  ####################################
  ### PART III: 3. Null IV GAM model
  ###              with NO exposure smooth
  ####################################
  if(messages == TRUE){ message("17. running full TSLS NULL GAM model") }
  iv_gam_mod0 = gamfit( wdata = mod_data,
                        dependent = outcome,
                        independent = NA,
                        linear_covariates = c(linear_covariates, "iv_predicted_exposure"),
                        smooth_covariates = smooth_covariates)

  ####################################
  ### PART III: 4. Full IV GAM model
  ###              with exposure smooth
  ####################################
  if(messages == TRUE){ message("16. running full TSLS GAM model") }
  iv_gam_mod = gamfit( wdata = mod_data,
                       dependent = outcome,
                       independent = "iv_predicted_exposure",
                       linear_covariates = linear_covariates,
                       smooth_covariates = smooth_covariates)

  ####################################
  ## PART III: 5. ANOVA F-test of
  ##              GAM_0 vs GAM
  ####################################
  if(messages == TRUE){ message("18. testing non-linearity of TSLS model: GAM vs NULL GAM") }
  a = anova(iv_gam_mod0$fit, iv_gam_mod$fit, test = "F")
  iv_nonlinearity_test = a[2,3:6]; names(iv_nonlinearity_test) = paste0( "", c("df","deviance","F","P"))

  nonlinearity_test = rbind(obs = obs_nonlinearity_test, tsls = iv_nonlinearity_test)

  ############################################
  ## PART III: 6. exposure summary statistics
  ############################################
  temp = na.omit( mod_data[, c(outcome, exposure, covariates, instrument)] )
  exp_ss = c( N = nrow(temp) ,
              mean = mean(temp[, exposure], na.rm = TRUE),
              min = min(temp[, exposure], na.rm = TRUE),
              max = max(temp[, exposure], na.rm = TRUE),
              sd = sd(temp[, exposure], na.rm = TRUE) )
  rm(temp)

  ####################################
  ### PART III: 7. Estimate instrument on exposure
  ###              FULL DATA Beta_ie
  ####################################
  if(messages == TRUE){ message("13. estimating effect of instrument on exposure for full data set") }
  beta_ie = iv_estimates( wdata = mod_data,
                         dependent = exposure,
                         instrument = instrument,
                         covariates = covariates)

  ####################################
  ### PART III: 8. Estimate instrument on outcome
  ###              FULL DATA Beta_io
  ####################################
  if(messages == TRUE){ message("13. estimating effect of instrument on exposure for full data set") }
  beta_io = iv_estimates( wdata = mod_data,
                          dependent = outcome,
                          instrument = instrument,
                          covariates = covariates)

  ####################################
  ### PART III: 9. Estimate causal estimate
  ###              FULL DATA Beta_iv
  ####################################
  if(messages == TRUE){ message("13. estimating effect of instrument on exposure for full data set") }
  ratio_mr_ss = as.data.frame( rbind( beta_ie = beta_ie,
                                      beta_io = beta_io,
                                      beta_iv = rep(NA, length(beta_io) ) ) )

  ratio_mr_ss$beta[3] = ratio_mr_ss$beta[2] / ratio_mr_ss$beta[1]
  ratio_mr_ss$se[3] = ratio_mr_ss$se[2] / ratio_mr_ss$beta[1]
  ratio_mr_ss$P[3] = ratio_mr_ss$P[2]
  ratio_mr_ss$n[3] = ratio_mr_ss$n[2]

  ratio_mr_ss = cbind(ratio_mr_ss, rbind(exp_ss, exp_ss, exp_ss))

  ####################################
  ### PART III: 10. ivreg()
  ###              FULL DATA Beta_iv
  ####################################
  if(messages == TRUE){ message("12. running a full TSLS model with ivreg") }
  iv_fit = ivregfit( wdata = mod_data,
                     outcome = outcome,
                     exposure = exposure,
                     instrument = instrument,
                     covariates = covariates,
                     weights_variable = NA,
                     rnt_outcome = FALSE)

  ####################################
  ### PART III: 11. Combine FULL DATA
  ###               TSLS | MR summary statistics
  ####################################
  ## ivreg sum stats
  ivreg_ss = c( iv_fit$summary, exp_ss )
  m = match(names(ratio_mr_ss), names(ivreg_ss) )
  tsls_ss = rbind( ratio_mr_ss, ivreg = ivreg_ss[m] )

  ####################################
  ## PART III: 12. Stratify data by
  ##               IV free exposure
  ####################################
  if(messages == TRUE){ message("9. stratifying the data by iv free exposure") }
  iv_strata_data = stratify_data( wdata = mod_data, stratify_on = "iv_free_exposure", strata = strata )

  ####################################
  ## PART III: 13. Summary Statistics
  ##               for exposure by
  ##               instrument free
  ##               exposure strata
  ####################################
  if(messages == TRUE){ message("10. stratified summary statistic for exposure") }
  iv_strata_exp_ss = stratify_sumstats(wdata = iv_strata_data, exposure = exposure)


  ####################################
  ### PART III: 14. Estimate instrument on exposure
  ###               beta coefficients by strata
  ####################################
  if(messages == TRUE){ message("13. estimating effect of instrument on exposure for each strata") }
  iv_strata_beta_ie = stratify_lmfit( wdata = iv_strata_data,
                                   outcome = exposure,
                                   exposure = instrument,
                                   covariates = covariates)
  iv_strata_beta_ie = cbind(iv_strata_beta_ie, iv_strata_exp_ss )

  ####################################
  ### PART III: 15. Estimate instrument on outcome
  ###               beta coefficients by strata
  ####################################
  if(messages == TRUE){ message("13. estimating effect of instrument on exposure for each strata") }
  iv_strata_beta_io = stratify_lmfit( wdata = iv_strata_data,
                                   outcome = outcome,
                                   exposure = instrument,
                                   covariates = covariates)
  iv_strata_beta_io = cbind(iv_strata_beta_io, iv_strata_exp_ss )

  ####################################
  ### PART III: 16. Beta_iv causal effect
  ###               estimates by strata
  ####################################
  iv_strata_ratio = stratify_ivratio( strata_beta_ie = iv_strata_beta_ie,
                    strata_beta_io = iv_strata_beta_io,
                    beta_ie = beta_ie["beta"],
                    strata_exp_ss = iv_strata_exp_ss,
                    tsls_ss = tsls_ss)

  ####################################
  ## PART III: 17. Run a ivreg model
  ##               on each strata
  ####################################
  if(messages == TRUE){ message("19. running ivreg on each strata") }
  iv_strata_ivreg = stratify_ivregfit( wdata = iv_strata_data,
                                outcome = outcome,
                                exposure = exposure,
                                instrument = instrument,
                                covariates = covariates,
                                weights_variable = NA,
                                rnt_outcome = FALSE)

  ####################################
  ## PART III: 18. Combine summary stats
  ##               for strafied ivreg and exposure
  ####################################
  iv_strata_ivreg = cbind(iv_strata_ivreg,  iv_strata_exp_ss)
  ## add the full data ivreg results
  iv_strata_ivreg = rbind(iv_strata_ivreg, fulldata = c( iv_fit$summary, exp_ss ))


  ####################################
  ## RESULTS OUT
  ####################################
  if(messages == TRUE){ message("22. compiling data to report") }
  names(rnt_outcome) = "outcome_RNTransformed"
  ## summary stats
  ss_out = data.frame( ## return input variable
                       outcome = outcome,
                       exposure = exposure,
                       instrument = instrument,
                       ## summary stats
                       number_of_outcome_outliers = number_of_outcome_outliers,
                       number_of_exposure_outliers = number_of_exposure_outliers,
                       W_outcome = W_outcome,
                       W_exposure = W_exposure,
                       rnt_outcome = rnt_outcome,
                       number_of_strata = length(strata_data),

                       stringsAsFactors = FALSE )
  rownames(ss_out) = "sumstats"

  ############################################
  ### Place all Summary Tables togethers
  ############################################
  summary_tables = list(
    ## test of non-linear relationships
    nonlinearity_test = nonlinearity_test,
    observational_sumstats = obs_ss,
    ### return the instrument exposure and instrument outcome stat estimates
    strata_beta_ie  = iv_strata_beta_ie,
    strata_beta_io  = iv_strata_beta_io,
    ## Complete data set TSLS estimates
    fulldata_tsls_sumstats = tsls_ss,
    ## Stratified TSLS estimates
    strata_tsls_ratio_sumstats = iv_strata_ratio,
    strata_tsls_ivreg_sumstats = iv_strata_ivreg,
    ## return the model data data frame
    model_data = mod_data
  )

  ############################################
  ### Place all complete data models in a single list
  ############################################
  models_no_stratification = list(
    ## observational models
    observational_gam_null = gam_mod0$fit,
    observational_gam = gam_mod$fit,
    observational_lm = lm_mod$fit,
    ## tsls models
    tsls_gam_null = iv_gam_mod0$fit,
    tsls_gam = iv_gam_mod$fit,
    tsls_ivreg = iv_fit$fit
  )

  ############################################
  ### Pull data together
  ############################################
  out = list(summary_stats = ss_out,
             summary_tables = summary_tables,
             models_no_stratification = models_no_stratification
             )

  if(messages == TRUE){ message("23. returning results to user") }
  return(out)

} ## end of function
