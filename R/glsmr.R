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

  ############################################
  ### I. Define Model Data
  ############################################
  # names(outcome) = "outcome"
  if(messages == TRUE){ message("1. Defining model data frame.") }

  model_variables = na.omit( c(outcome, exposure, instrument, linear_covariates, smooth_covariates, weights_variable) )
  mod_data = wdata[, c(model_variables)]

  ############################################
  ## II. Identify outcome outliers
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
  ## III. Identify exposure outliers
  ############################################
  outliers = id_outliers( y = mod_data[, exposure], outlier_method = outlier_method, outlier_cutoff = outlier_cutoff)

  # How many exposure outliers
  number_of_exposure_outliers = length(outliers); names(number_of_exposure_outliers) = "number_of_exposure_outliers"

  ## Turn outliers to NA
  if(length(outliers) > 0){
    mod_data[outliers, exposure] = NA
  }

  ############################################
  ## IV. normality of outcome
  ############################################
  if(messages == TRUE){ message("3. Estimate Shapiro-Wilk normality W-stat") }
  W_outcome = normW(mod_data[, outcome]); names(W_outcome) = "W_outcome"

  ############################################
  ## IV. normality of exposure
  ############################################
  W_exposure = normW(mod_data[, exposure]); names(W_exposure) = "W_exposure"

  ############################################
  ## V. rank normalize the outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    if(messages == TRUE){ message("4. Performing rank normal transformation of outcome") }
    mod_data[, outcome] = rntransform( mod_data[, outcome] )
  } else {
    if(messages == TRUE){ message("4. No rank normal transformation of outcome performed") }
  }

  ######################
  ### VI. Linear model
  ######################
  if(messages == TRUE){ message("5. running full observational linear model") }
  cvs = na.omit(c(linear_covariates, smooth_covariates))
  if(length(cvs) == 0){cvs = NA}
  lm_mod = lmfit( wdata = mod_data,
                  dependent = outcome,
                  independent = exposure,
                  covariates = cvs)

  ######################
  ### VII. GAM model
  ######################
  if(messages == TRUE){ message("6. running full observational GAM model") }
  gam_mod = gamfit( wdata = mod_data,
                  dependent = outcome,
                  independent = exposure,
                  linear_covariates = linear_covariates,
                  smooth_covariates = smooth_covariates)

  ####################################
  ### VII. GAM model with no exposure smooth
  ###      exposure is just another parametric
  ###      parameter.
  ####################################
  if(messages == TRUE){ message("7. running full observational NULL GAM model") }
  gam_mod0 = gamfit( wdata = mod_data,
                    dependent = outcome,
                    independent = NA,
                    linear_covariates = c(linear_covariates, exposure),
                    smooth_covariates = smooth_covariates)


  ####################################
  ## VIII. ANOVA test of GAM with no
  ##       exposure smooth and full GAM
  ####################################
  if(messages == TRUE){ message("8. testing non-linearity of observational data: GAM vs NULL GAM") }
  a = anova(gam_mod0$fit, gam_mod$fit, test = "F")
  obs_nonlinearity_test = a[2,3:6]; names(obs_nonlinearity_test) = paste0( "obs_nonlinearity_test_", c("df","deviance","F","P"))

  ####################################
  ## IX. Stratify
  ####################################
  if(messages == TRUE){ message("9. stratifying the data") }
  if(length(strata) == 1){
    q = quantile( mod_data[, exposure], probs = seq(0, 1, 1/strata), na.rm = TRUE )
  } else {
    if(class(strata) == "numeric" & length(strata) >=3 ){
      q = strata
    } else {
      stop("please check strata parameter. Acceptable values are a single numeric value indicating the number of quantiles, or a numeric vector of at least length 3 to define strata boundries.")
    }
  }

  ### number of strata
  number_of_strata = length(q)-1; names(number_of_strata) = "number_of_strata"

  ## identify samples of each strata
  ## and add a strata label to the model data frame
  mod_data$strata = as.factor( cut(mod_data[, exposure], q, include.lowest=T, labels=F) )

  ### make a data frame for each strata
  strata_data = lapply(1:number_of_strata, function(i){
    w = which(mod_data$strata == i)
    mod_data[w,]
  })


  ## estimate exposure summary statistics for samples in each strata
  if(messages == TRUE){ message("10. estimating strata summary statistics (n, mean, min, max).") }
  exposure_strata_sumstats = t( sapply( 1:length(strata_data), function(i){
    x = strata_data[[i]]
    N = length(x[, exposure])
    r = range(x[, exposure], na.rm = TRUE)
    m = mean(x[, exposure], na.rm = TRUE)
    s = sd(x[, exposure], na.rm = TRUE)
    out = c(N, m, r, s); names(out) = c("N","mean","min","max", "sd")
    return(out)
  }) )

  ####################################
  ## X. Run a linear model on each
  ##    strata
  ####################################
  if(messages == TRUE){ message("11. running observational linear models on each strata") }
  cvs = na.omit(c(linear_covariates, smooth_covariates))

  strata_linear_mods = t( sapply(1:length(strata_data), function(i){
    x = strata_data[[i]]
    ### variability check
    var_present = apply(x, 2, function(i){ length(unique(i))}) ## only column name 'strata' should have a var == 1
    if( sum(var_present == 1) == 1){
      lm_mod = lmfit( wdata = x,
                           dependent = outcome,
                           independent = exposure,
                           covariates = cvs)
    }else{
      o = rep(NA, 10); names(o) = c("n","W","Rsq","Fstat","df","dendf","beta","se","tval","P")
      lm_mod = list(summary = o)
    }

    return(lm_mod$summary)
  }) )

  ####################################
  ## XI. Add strata sum stats to
  ##    linear model estimates
  ####################################
  strata_linear_mods = cbind(strata_linear_mods, exposure_strata_sumstats)
  rownames(strata_linear_mods) = paste0("strata_", 1:nrow(strata_linear_mods))
  strata_linear_mods = as.data.frame(strata_linear_mods)


  ######################
  ### XII. IVreg model
  ######################
  if(messages == TRUE){ message("12. running a full TSLS model with ivreg") }
  cvs = na.omit(c(linear_covariates, smooth_covariates))

  iv_fit = ivregfit( wdata = mod_data,
            outcome = outcome,
            exposure = exposure,
            instrument = instrument,
            covariates = cvs,
            weights_variable = NA,
            rnt_outcome = FALSE)

  ######################
  ### XIV. Estimate d.hat
  ###      AND beta coefficents
  ######################
  if(messages == TRUE){ message("13. estimating d.hat (a.k.a IV predicted exposure)") }
  ## Univariate analysis
  form = formula(paste0(exposure, " ~ ", instrument ))
  fit_grs_exp = lm( form , data = mod_data)

  ## Variance in Exposure Explained by Instrument
  Exposure_Var_Exp_by_IV = summary(fit_grs_exp)$r.sq; names(Exposure_Var_Exp_by_IV) = "Exposure_Var_Exp_by_IV"

  ## d.hat
  temp = fitted.values(fit_grs_exp)
  m = match(rownames(mod_data), names(temp))
  mod_data$d.hat = temp[m]

  ## Univariate Coefficients
  if(messages == TRUE){ message("14. estimating univariate variance explained for IV on exposure") }
  univariate = summary(fit_grs_exp)$coefficients[instrument,c(1,2,4)]
  names( univariate ) = paste0( "grs_on_exp_", c("beta","se","P") )

  ## Multivariate estimate
  if(messages == TRUE){ message("15. estimating multivariate variance explained for IV on exposure") }
  cvs = na.omit( c(linear_covariates, smooth_covariates) )
  if( length( cvs ) > 0 ) {
    form = formula(paste0(exposure, " ~ ", paste0( cvs, collapse = " + ") , " + ", instrument ))
    fit_grs_exp = lm( form , data = mod_data)
    multivariate = summary(fit_grs_exp)$coefficients[instrument,c(1,2,4)]
    ## variance explained given the multivariate model
    a = anova(fit_grs_exp)
    etasq = a[,2]/sum(a[,2]); names(etasq) = paste0("etasq_",rownames(a))
    Exposure_Var_Exp_by_IV_etasq = etasq[ paste0("etasq_", instrument )]; names(Exposure_Var_Exp_by_IV_etasq) = "Exposure_Var_Exp_by_IV_etasq"
    } else {
      multivariate = c(NA, NA, NA)
      Exposure_Var_Exp_by_IV_etasq = NA; names(Exposure_Var_Exp_by_IV_etasq) = "Exposure_Var_Exp_by_IV_etasq"
    }
  names( multivariate ) = paste0( "grs_on_exp_", c("beta","se","P") )

  ## combine univariate and multivariate estimates
  iv_on_exp_coeff = rbind(univariate, multivariate)


  ######################
  ### XIV. IV GAM model
  ######################
  if(messages == TRUE){ message("16. running full TSLS GAM model") }
  iv_gam_mod = gamfit( wdata = mod_data,
                    dependent = outcome,
                    independent = "d.hat",
                    linear_covariates = linear_covariates,
                    smooth_covariates = smooth_covariates)

  ####################################
  ### XV. IV GAM model with no exposure smooth
  ###      exposure is just another parametric
  ###      parameter.
  ####################################
  if(messages == TRUE){ message("17. running full TSLS NULL GAM model") }
  iv_gam_mod0 = gamfit( wdata = mod_data,
                     dependent = outcome,
                     independent = NA,
                     linear_covariates = c(linear_covariates, "d.hat"),
                     smooth_covariates = smooth_covariates)


  ####################################
  ## XVI ANOVA test of GAM with no
  ##       exposure smooth and full GAM
  ####################################
  if(messages == TRUE){ message("18. testing non-linearity of TSLS model: GAM vs NULL GAM") }
  a = anova(iv_gam_mod0$fit, iv_gam_mod$fit, test = "F")
  iv_nonlinearity_test = a[2,3:6]; names(iv_nonlinearity_test) = paste0( "iv_nonlinearity_test_", c("df","deviance","F","P"))

  ####################################
  ## XVII. Run a ivreg model on each
  ##      strata
  ####################################
  if(messages == TRUE){ message("19. running ivreg on each strata") }
  cvs = na.omit(c(linear_covariates, smooth_covariates))

  strata_ivreg_mods = t( sapply(1:length(strata_data), function(i){
    x = strata_data[[i]]
    ### variability check
    var_present = apply(x, 2, function(i){ length(unique(i))})
    if( sum(var_present == 1) == 1){

      iv_mod = ivregfit( wdata = x,
                         outcome = outcome,
                         exposure = exposure,
                         instrument = instrument,
                         covariates = cvs,
                         weights_variable = NA,
                         rnt_outcome = FALSE)
    }else{
      o = rep(NA, 18); names(o) = c("n","Rsq","Wald_stat","Wald_df1","Wald_df2","Wald_P",
                                    "Fstat","F_df1","F_df2","F_P","WuH_stat","WuH_df1","WuH_df2","WuH_P",
                                    "beta","se","tval","P")
      iv_mod = list(summary = o)
    }

    return(iv_mod$summary)
  }) )

  ####################################
  ## XVIII. Add strata sum stats to
  ##    linear model estimates
  ####################################
  strata_ivreg_mods = cbind(strata_ivreg_mods, exposure_strata_sumstats)
  rownames(strata_ivreg_mods) = paste0("strata_", 1:nrow(strata_ivreg_mods))
  strata_ivreg_mods = as.data.frame(strata_ivreg_mods)

  ####################################
  ## XIX Run a linear model of
  ##      instrument on outcome
  ##       for each strata
  ####################################
  if(messages == TRUE){ message("20. running a linear model of IV on outcome for each strata") }
  cvs = na.omit(c(linear_covariates, smooth_covariates))

  strata_IV_linear_mods = t( sapply(1:length(strata_data), function(i){
    x = strata_data[[i]]
    ### variability check
    var_present = apply(x, 2, function(i){ length(unique(i))}) ## only column name 'strata' should have a var == 1
    if( sum(var_present == 1) == 1){
      lm_mod = lmfit( wdata = x,
                      dependent = outcome,
                      independent = instrument,
                      covariates = cvs)
    }else{
      o = rep(NA, 10); names(o) = c("n","W","Rsq","Fstat","df","dendf","beta","se","tval","P")
      lm_mod = list(summary = o)
    }

    return(lm_mod$summary)
  }) )
  strata_IV_linear_mods = as.data.frame(strata_IV_linear_mods)

  ## derive ratio values
  if(messages == TRUE){ message("21. deriving MR ratio for each strata") }
  if( !is.na( iv_on_exp_coeff["multivariate", "grs_on_exp_beta"] ) ){
    strata_IV_linear_mods$beta_ratio = strata_IV_linear_mods$beta / iv_on_exp_coeff["multivariate", "grs_on_exp_beta"]
    strata_IV_linear_mods$se_ratio = strata_IV_linear_mods$se / iv_on_exp_coeff["multivariate", "grs_on_exp_beta"]
  } else {
    strata_IV_linear_mods$beta_ratio = strata_IV_linear_mods$beta / iv_on_exp_coeff["univariate", "grs_on_exp_beta"]
    strata_IV_linear_mods$se_ratio = strata_IV_linear_mods$se / iv_on_exp_coeff["univariate", "grs_on_exp_beta"]
  }

  ####################################
  ## XI. Add strata sum stats to
  ##    linear model estimates
  ####################################
  strata_IV_linear_mods = cbind(strata_IV_linear_mods, exposure_strata_sumstats)
  rownames(strata_IV_linear_mods) = paste0("strata_", 1:nrow(strata_IV_linear_mods))
  strata_IV_linear_mods = as.data.frame(strata_IV_linear_mods)

  ####################################
  ## RESULTS OUT
  ####################################
  if(messages == TRUE){ message("22. compiling data to report") }
  names(rnt_outcome) = "outcome_RNTransformed"
  ## summary stats
  ss_out = data.frame( number_of_outcome_outliers = number_of_outcome_outliers,
                       number_of_exposure_outliers = number_of_exposure_outliers,
                       W_outcome = W_outcome,
                       W_exposure = W_exposure,
                       rnt_outcome = rnt_outcome,
                       number_of_strata = number_of_strata,
                       obs_nonlinearity_test_df = obs_nonlinearity_test[1],
                       obs_nonlinearity_test_deviance = obs_nonlinearity_test[2],
                       obs_nonlinearity_test_F = obs_nonlinearity_test[3],
                       obs_nonlinearity_test_P = obs_nonlinearity_test[4],
                       exposure_VarExp_by_iv = Exposure_Var_Exp_by_IV,
                       exposure_VarExp_by_iv_etasq = Exposure_Var_Exp_by_IV_etasq,
                       iv_nonlinearity_test_df = iv_nonlinearity_test[1],
                       iv_nonlinearity_test_deviance = iv_nonlinearity_test[2],
                       iv_nonlinearity_test_F = iv_nonlinearity_test[3],
                       iv_nonlinearity_test_P = iv_nonlinearity_test[4],
                       exposure = exposure,
                       outcome = outcome, 
                       stringsAsFactors = FALSE )
  rownames(ss_out) = "sumstats"

  out = list(strata_linear_mods = strata_linear_mods,
             strata_ivreg_mods = strata_ivreg_mods,
             strata_IV_linear_mods = strata_IV_linear_mods,
             summary_stats = ss_out,
             full_linear_model = lm_mod[[1]],
             null_full_gam_model = gam_mod0[[1]],
             full_gam_model = gam_mod[[1]],
             iv_on_exp_coeff = iv_on_exp_coeff,
             full_ivreg_model = iv_fit[[1]],
             null_full_iv_gam_model = iv_gam_mod0[[1]],
             full_iv_gam_model = iv_gam_mod[[1]],
             model_data = mod_data
             )
  if(messages == TRUE){ message("23. returning results to user") }
  return(out)

} ## end of function
