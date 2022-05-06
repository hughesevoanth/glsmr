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
#' @param strata a single integer or vector to define how strata should be defined.
#' * `strata = 4`: will define quartiles or 4 evenly sized bins
#' * `strata = 10`: will define deciles or 10 evenly sized bins
#' * `strata = c(1,10,20,30)`: a user defined numeric vector to define boundaries for each strata.
#'    - The numeric vector example provided will define 4 strata. Lower bound values are inclusive, upper bounds are exclusive, to the exception of the last bounding value.
#' @param rnt_outcome binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
#' @param weights_variable a single string character of the column name for a weights variable
#' @param outlier_cutoff a single numeric value to define a cutoff value for how many iqr or sd units outlier values
#' @param outlier_method a single string character of "iqr" or "sd" to determine if outlier should be determined by means and sd or medians and iqr.
#' @param messages should a progress message be printed to screen - binary TRUE or FALSE
#' @param return_models should the model data data frame and each gam and linear model be returned (TRUE or FALSE). Default is FALSE.
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
                  messages = FALSE,
                  return_models = FALSE){

  ############################################
  ### 0. look for any errors in parameters
  ############################################
  # if (!strata[1] %in% c("quartiles","deciles") & length(strata) < 3 & class(strata) != "numeric" ){
  #   stop("strata parameter can either be defined as (1) quartiles, (2) deciles, or (3) a numeric vector of at least length 3 defining strata boundries")
  #  }

  ## PART I: SETTING UP THE DATA
  if(messages == TRUE){ message("Part I. Setting up the data.") }

  ############################################
  ### PART I: 1. Define Model Data
  ############################################
  if(messages == TRUE){ message("Part I.1. Defining model data frame.") }
  model_variables = na.omit( c(outcome, exposure, instrument, linear_covariates, smooth_covariates, weights_variable) )
  # mod_data = na.omit( wdata[, c(model_variables)] )
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

  ############################################
  ## PART I: 8. Derive IV-Free Exposure
  ############################################
  if(messages == TRUE){ message("Part I.8. derive instument free exposure") }
  mod_data = iv_free_exposure( wdata = mod_data,
                               exposure = exposure,
                               instrument = instrument,
                               covariates = NA,
                               exposure_mean_normalize = TRUE)

  ####################################
  ## PART I: 9. Derive Instrument predicted exposure or d.hat
  ####################################
  if(messages == TRUE){ message("Part I.9. derive instument predicted exposure") }
  mod_data = iv_predicted_exposure( wdata = mod_data,
                                    exposure = exposure,
                                    instrument = instrument,
                                    covariates = covariates)

  ############################################
  ## PART I: 10. Summary Statistics for exposure
  ############################################
  if(messages == TRUE){ message("Part I.10. summary statistic for exposure") }
  temp = na.omit( mod_data[, c(outcome, exposure, covariates)] )
  exp_ss = c( N = nrow(temp) ,
              mean = mean(temp[, exposure], na.rm = TRUE),
              min = min(temp[, exposure], na.rm = TRUE),
              max = max(temp[, exposure], na.rm = TRUE),
              sd = sd(temp[, exposure], na.rm = TRUE) )
  rm(temp)

  ### Estimate again with filtering on instruemtn NAs too
  temp = na.omit( mod_data[, c(outcome, exposure, covariates, instrument)] )
  iv_exp_ss = c( N = nrow(temp) ,
              mean = mean(temp[, exposure], na.rm = TRUE),
              min = min(temp[, exposure], na.rm = TRUE),
              max = max(temp[, exposure], na.rm = TRUE),
              sd = sd(temp[, exposure], na.rm = TRUE) )
  rm(temp)

  ############################################
  ## PART I: 11. Stratify data by Exposure
  ############################################
  if(messages == TRUE){ message("Part I.11. stratifying the data by the iv free exposure") }
  strata_data = stratify_data( wdata = mod_data, stratify_on = "iv_free_exposure", strata = strata )

  ############################################
  ## PART I: 12. Summary Statistics for exposure
  ##             by strata
  ############################################
  if(messages == TRUE){ message("Part I.12. stratified summary statistic for exposure") }
  strata_exp_ss = stratify_sumstats(wdata = strata_data, exposure = exposure, model_vars = c(outcome, exposure, covariates) )
  iv_strata_exp_ss = stratify_sumstats(wdata = strata_data, exposure = exposure, model_vars = c(outcome, exposure, covariates, instrument) )

  ## PART II: OBSERVATIONAL MODELING
  if(messages == TRUE){ message("Part II. Observational modeling") }
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

  ############################################
  ### PART II: 2. full GAM model
  ############################################
  if(messages == TRUE){ message("Part II.2. running full observational GAM model") }
  gam_mod = gamfit( wdata = mod_data,
                    dependent = outcome,
                    independent = exposure,
                    linear_covariates = linear_covariates,
                    smooth_covariates = smooth_covariates)


  ############################################
  ## PART II: 3. ANOVA test of GAM with no
  ##             exposure smooth and full GAM
  ############################################
  if(messages == TRUE){ message("Part II.3. testing non-linearity of observational data: GAM vs NULL GAM") }
  a = anova(gam_mod0$fit, gam_mod$fit, test = "F")
  obs_nonlinearity_test = a[2,3:6]; names(obs_nonlinearity_test) = paste0( "", c("df","deviance","F","P"))


  ############################################
  ### PART II: 4. Linear model
  ############################################
  if(messages == TRUE){ message("Part II.4. running full observational linear model") }
  lm_mod = lmfit( wdata = mod_data,
                  dependent = outcome,
                  independent = exposure,
                  covariates = covariates)

  ## Combine summary stats of full for data observational data
  obs_ss = c(lm_mod$summary, exp_ss)

  ############################################
  ## PART II: 5. Run observational linear model on each
  ##             strata
  ############################################
  if(messages == TRUE){ message("Part II.5. running observational linear models on each strata") }
  strata_lm_mod = stratify_lmfit( wdata = strata_data,
                               outcome = outcome,
                               exposure = exposure,
                               covariates = covariates)

  ## Combine strata sumstats and linear model estimates
  strata_obs_ss = cbind(strata_lm_mod, strata_exp_ss)

  ############################################
  ## PART II: 6. Observational Meta
  ##              Analysis
  ############################################
  if(messages == TRUE){ message("Part II.6. run meta analysis on observational strata betas and means") }
  obs_meta_test = meta_test( strata_obs_ss )

  ############################################
  ## PART II: 7. Combine strata summary stats
  ##              & full model summary stats
  ############################################
  obs_ss = rbind( strata_obs_ss , fulldata = obs_ss )
  obs_ss = obs_ss[, -c(4:6,9) ]

  ## PART III: OBSERVATIONAL MODELING
  if(messages == TRUE){ message("Part III. MR modeling") }
  ####################################
  ### PART III: 1. Null IV GAM model
  ###              with NO exposure smooth
  ####################################
  if(messages == TRUE){ message("Part III.1. running full TSLS NULL GAM model") }
  iv_gam_mod0 = gamfit( wdata = mod_data,
                        dependent = outcome,
                        independent = NA,
                        linear_covariates = c(linear_covariates, "iv_predicted_exposure"),
                        smooth_covariates = smooth_covariates)

  ####################################
  ### PART III: 4. Full IV GAM model
  ###              with exposure smooth
  ####################################
  if(messages == TRUE){ message("Part III.2. running full TSLS GAM model") }
  iv_gam_mod = gamfit( wdata = mod_data,
                       dependent = outcome,
                       independent = "iv_predicted_exposure",
                       linear_covariates = linear_covariates,
                       smooth_covariates = smooth_covariates)

  ####################################
  ## PART III: 3. ANOVA F-test of
  ##              GAM_0 vs GAM
  ####################################
  if(messages == TRUE){ message("Part III.3. testing non-linearity of TSLS model: GAM vs NULL GAM") }
  a = anova(iv_gam_mod0$fit, iv_gam_mod$fit, test = "F")
  iv_nonlinearity_test = a[2,3:6]; names(iv_nonlinearity_test) = paste0( "", c("df","deviance","F","P"))

  nonlinearity_test = rbind(obs = obs_nonlinearity_test, tsls = iv_nonlinearity_test)

  ####################################
  ### PART III: 4. Estimate instrument on exposure
  ###              FULL DATA Beta_ie
  ####################################
  if(messages == TRUE){ message("Part III.4. estimating effect of instrument on exposure (beta_ie) for full data set") }
  beta_ie = iv_estimates( wdata = mod_data,
                         dependent = exposure,
                         instrument = instrument,
                         covariates = covariates)

  ####################################
  ### PART III: 5. Estimate instrument on outcome
  ###              FULL DATA Beta_io
  ####################################
  if(messages == TRUE){ message("Part III.5. estimating effect of instrument on outcome (beta_io) for full data set") }
  beta_io = iv_estimates( wdata = mod_data,
                          dependent = outcome,
                          instrument = instrument,
                          covariates = covariates)

  ####################################
  ### PART III: 6. Estimate causal estimate
  ###              FULL DATA Beta_iv
  ####################################
  if(messages == TRUE){ message("Part III.6. estimating causal effect estimate (beta_iv)") }
  ratio_mr_ss = as.data.frame( rbind( beta_ie = beta_ie,
                                      beta_io = beta_io,
                                      beta_iv = rep(NA, length(beta_io) ) ) )

  ratio_mr_ss$beta[3] = ratio_mr_ss$beta[2] / ratio_mr_ss$beta[1]
  ratio_mr_ss$se[3] = ratio_mr_ss$se[2] / ratio_mr_ss$beta[1]
  ratio_mr_ss$P[3] = ratio_mr_ss$P[2]
  ratio_mr_ss$n[3] = ratio_mr_ss$n[2]

  ratio_mr_ss = cbind(ratio_mr_ss, rbind(iv_exp_ss, iv_exp_ss, iv_exp_ss))

  ####################################
  ### PART III: 7. ivreg()
  ###              FULL DATA Beta_iv
  ####################################
  if(messages == TRUE){ message("Part III.7. running a full TSLS model with ivreg") }
  iv_fit = ivregfit( wdata = mod_data,
                     outcome = outcome,
                     exposure = exposure,
                     instrument = instrument,
                     covariates = covariates,
                     weights_variable = NA,
                     rnt_outcome = FALSE)

  ### Combine FULL DATA TSLS | MR summary statistics
  ivreg_ss = c( iv_fit$summary, iv_exp_ss )
  m = match(names(ratio_mr_ss), names(ivreg_ss) )
  tsls_ss = rbind( ratio_mr_ss, ivreg = ivreg_ss[m] )


  ####################################
  ### PART III: 8. Estimate instrument on exposure
  ###               beta coefficients by strata
  ####################################
  if(messages == TRUE){ message("Part III.8. estimating effect of instrument on exposure for each strata") }
  iv_strata_beta_ie = stratify_lmfit( wdata = strata_data,
                                   outcome = exposure,
                                   exposure = instrument,
                                   covariates = covariates)
  iv_strata_beta_ie = cbind(iv_strata_beta_ie, iv_strata_exp_ss )

  ####################################
  ### PART III: 9. Testing assumption
  ###                that IV-exposure effect is
  ###                stable across strata
  ###                i.e. relationship is linear
  ####################################
  if(messages == TRUE){ message("Part III.9. IV-exposure assumption meta (beta_ie)") }
  ie_meta_test = meta_test( iv_strata_beta_ie )


  ####################################
  ### PART III: 10. Estimate instrument on outcome
  ###               beta coefficients by strata
  ####################################
  if(messages == TRUE){ message("Part III.10. estimating effect of instrument on outcome for each strata") }
  iv_strata_beta_io = stratify_lmfit( wdata = strata_data,
                                   outcome = outcome,
                                   exposure = instrument,
                                   covariates = covariates)
  iv_strata_beta_io = cbind(iv_strata_beta_io, iv_strata_exp_ss )

  ####################################
  ### PART III: 11. Beta_iv causal effect
  ###               estimates by strata
  ####################################
  if(messages == TRUE){ message("Part III.11. estimating the iv ratio MR estimates for each strata") }
  iv_strata_ratio = stratify_ivratio( strata_beta_ie = iv_strata_beta_ie,
                    strata_beta_io = iv_strata_beta_io,
                    beta_ie = beta_ie["beta"],
                    strata_exp_ss = iv_strata_exp_ss,
                    tsls_ss = tsls_ss)

  ####################################
  ## PART III: 12. Run a ivreg model
  ##               on each strata
  ####################################
  if(messages == TRUE){ message("Part III.12. running ivreg on each strata") }
  iv_strata_ivreg = stratify_ivregfit( wdata = strata_data,
                                outcome = outcome,
                                exposure = exposure,
                                instrument = instrument,
                                covariates = covariates,
                                weights_variable = NA,
                                rnt_outcome = FALSE)
  iv_strata_ivreg = cbind(iv_strata_ivreg, iv_strata_exp_ss )

  ####################################
  ### PART III: 13. ivreg meta
  ####################################
  if(messages == TRUE){ message("Part III.13. ivreg meta analysis (beta_iv)") }
  iv_meta_test = meta_test( iv_strata_ivreg )

  ## add the full data ivreg results
  iv_strata_ivreg = rbind(iv_strata_ivreg, fulldata = c( iv_fit$summary, iv_exp_ss ))


  if(messages == TRUE){ message("Part IV. returning results to user") }
  ####################################
  ## RESULTS OUT
  ####################################
  if(messages == TRUE){ message("Part IV.1. compiling data to report") }
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

  ## meta analysis stats.
  meta_stats = rbind(obs_meta_test, ie_meta_test, iv_meta_test)
  rownames(meta_stats) = c("obs","ie","mr")

  ############################################
  ### Place all Summary Tables together
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
    ## Stratified meta sum stats
    strata_meta_analysis = meta_stats,
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
  if(return_models == TRUE){
  out = list(summary_stats = ss_out,
             summary_tables = summary_tables,
             models_no_stratification = models_no_stratification
             )
  } else {
    out = list(summary_stats = ss_out,
               summary_tables = summary_tables
    )
  }

  if(messages == TRUE){ message("Part IV.2. returning results to user") }
  return(out)

} ## end of function
