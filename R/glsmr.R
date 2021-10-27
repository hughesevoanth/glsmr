#' GAM and linear-strata Mendelian randomization
#'
#' This function derives strata defined observational and TSLS (MR) estimates using linear modeling.
#' The full data set will be fit to a linear model and a generalized linear model (GAM).
#' Then subsequently the exposure data will be stratified and a linear model will be used to derive
#' effect estimmates for each strata.
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
#' @param sd_outlier_cutoff a single numeric value to define a cutoff value for how many iqr or sd units outlier values
#' @param sd_or_iqr_cutoff a single string character of "iqr" or "sd" to determine if outlier should be determined by means and sd or medians and iqr.
#' @return returns a glsmr object containing the complete linear and GAM models for the full data set, summary statistics for the data, a strata observational table, and a strata TSLS (MR) table.
#' @importFrom stats na.omit anova quantile sd formula lm fitted.values
#' @export
#' @examples
#' glsmr()
glsmr = function( wdata,
                  outcome = NA,
                  exposure = NA,
                  instrument = NA,
                  linear_covariates = NA,
                  smooth_covariates = NA,
                  strata = "quartiles",
                  rnt_outcome = FALSE,
                  weights_variable = NA,
                  sd_outlier_cutoff = 5,
                  sd_or_iqr_cutoff = "iqr"){

  ############################################
  ### 0. look for any errors in parameters
  ############################################
  if (!strata[1] %in% c("quartiles","deciles") & length(strata) < 3 & class(strata) != "numeric" ){
    stop("strata parameter can either be defined as (1) quartiles, (2) deciles, or (3) a numeric vector of at least length 3 defining strata boundries")
    }

  ############################################
  ### I. Define Model Data
  ############################################
  # names(outcome) = "outcome"
  ##
  model_variables = na.omit( c(outcome, exposure, instrument, linear_covariates, smooth_covariates, weights_variable) )
  mod_data = wdata[, c(model_variables)]

  ############################################
  ## II. Identify outcome outliers
  ############################################
  outliers = id_outliers( y = mod_data[, outcome], iqr_or_sd = sd_or_iqr_cutoff, number_of_sd_iqr = sd_outlier_cutoff)

  # How many outcome outliers
  number_of_outcome_outliers = length(outliers); names(number_of_outcome_outliers) = "number_of_outcome_outliers"

  ## Turn outliers to NA
  if(length(outliers) > 0){
    mod_data[outliers, outcome] = NA
  }

  ############################################
  ## III. Identify exposure outliers
  ############################################
  outliers = id_outliers( y = mod_data[, exposure], iqr_or_sd = sd_or_iqr_cutoff, number_of_sd_iqr = sd_outlier_cutoff)

  # How many exposure outliers
  number_of_exposure_outliers = length(outliers); names(number_of_exposure_outliers) = "number_of_exposure_outliers"

  ## Turn outliers to NA
  if(length(outliers) > 0){
    mod_data[outliers, exposure] = NA
  }

  ############################################
  ## IV. normality of outcome
  ############################################
  W_outcome = normW(mod_data[, outcome]); names(W_outcome) = "W_outcome"

  ############################################
  ## IV. normality of exposure
  ############################################
  W_exposure = normW(mod_data[, exposure]); names(W_exposure) = "W_exposure"

  ############################################
  ## V. rank normalize the outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    mod_data[, outcome] = rntransform( mod_data[, outcome] )
  }

  ######################
  ### VI. Linear model
  ######################
  cvs = na.omit(c(linear_covariates, smooth_covariates))
  lm_mod = lmfit( wdata = mod_data,
                  dependent = outcome,
                  independent = exposure,
                  covariates = cvs)

  ######################
  ### VII. GAM model
  ######################
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
  gam_mod0 = gamfit( wdata = mod_data,
                    dependent = outcome,
                    independent = NA,
                    linear_covariates = c(linear_covariates, exposure),
                    smooth_covariates = smooth_covariates)


  ####################################
  ## VIII. ANOVA test of GAM with no
  ##       exposure smooth and full GAM
  ####################################
  a = anova(gam_mod0$fit, gam_mod$fit, test = "F")
  linear_nonlinear_test = a[2,3:6]; names(linear_nonlinear_test) = paste0( "linear_nonlinear_test_", c("df","deviance","F","P"))

  ####################################
  ## IX. Stratify
  ####################################
  if(strata[1] == "quartiles"){
    q = quantile( mod_data[, exposure], na.rm = TRUE )
  } else {
    if(strata[1] == "deciles"){
      q = quantile( mod_data[, exposure], probs = seq(0,1, by=0.1), na.rm = TRUE )
    } else {
      if(class(strata) == "numeric" & length(strata) >=3 ){
          q = strata
      } else {
        stop("please check strata parameter. Acceptable values are 'quartiles', 'deciles', or a numeric vector of at least length 3 to define strata boundries.")
      }
    }
  }

  ### identify samples of each strata
  number_of_strata = length(q)-1; names(number_of_strata) = "number_of_strata"
  strata_samples = lapply(1:number_of_strata, function(i){
    if(i != number_of_strata){
      w = which( mod_data[, exposure] >= q[i] & mod_data[, exposure] < q[i+1])
    } else {
      w = which( mod_data[, exposure] >= q[i] & mod_data[, exposure] <= q[i+1])
    }
    return(w)
  })
  ### make a data frame for each strata
  strata_data = lapply(strata_samples, function(i){
    mod_data[i,]
  })

  ## add a label to the model data frame
  mod_data$strata = 0
  for(i in 1:length(strata_samples)){
    w = strata_samples[[i]]
    mod_data$strata[w] = paste0("strata_", i)
  }
  mod_data$strata = as.factor(mod_data$strata)

  ## estimate exposure summary statistics for samples in each strata
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
  cvs = na.omit(c(linear_covariates, smooth_covariates))

  strata_linear_mods = t( sapply(1:length(strata_data), function(i){
    x = strata_data[[i]]
    ### variability check
    var_present = apply(x, 2, function(i){ length(unique(i))})
    if( sum(var_present == 1) == 0){
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
  ######################
  form = formula(paste0(exposure , "~", instrument))
  fitDhat = lm( form ,data = mod_data)
  VarExp_by_Instrument_on_Exposure = summary(fitDhat)$r.sq; names(VarExp_by_Instrument_on_Exposure) = "VarExp_by_Instrument_on_Exposure"
  temp = fitted.values(fitDhat)
  m = match(rownames(mod_data), names(temp))
  mod_data$d.hat = temp[m]

  ######################
  ### XIV. IV GAM model
  ######################
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
  iv_gam_mod0 = gamfit( wdata = mod_data,
                     dependent = outcome,
                     independent = NA,
                     linear_covariates = c(linear_covariates, "d.hat"),
                     smooth_covariates = smooth_covariates)


  ####################################
  ## XVI ANOVA test of GAM with no
  ##       exposure smooth and full GAM
  ####################################
  a = anova(iv_gam_mod0$fit, iv_gam_mod$fit, test = "F")
  iv_linear_nonlinear_test = a[2,3:6]; names(iv_linear_nonlinear_test) = paste0( "iv_linear_nonlinear_test_", c("df","deviance","F","P"))

  ####################################
  ## XVII. Run a ivreg model on each
  ##      strata
  ####################################
  cvs = na.omit(c(linear_covariates, smooth_covariates))

  strata_ivreg_mods = t( sapply(1:length(strata_data), function(i){
    x = strata_data[[i]]
    ### variability check
    var_present = apply(x, 2, function(i){ length(unique(i))})
    if( sum(var_present == 1) == 0){

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
  ## XI. Add strata sum stats to
  ##    linear model estimates
  ####################################
  strata_ivreg_mods = cbind(strata_ivreg_mods, exposure_strata_sumstats)
  rownames(strata_ivreg_mods) = paste0("strata_", 1:nrow(strata_ivreg_mods))
  strata_ivreg_mods = as.data.frame(strata_ivreg_mods)

  ####################################
  ## RESULTS OUT
  ####################################
  names(rnt_outcome) = "outcome_RNTransformed"
  ## summary stats
  ss_out = unlist( c(number_of_outcome_outliers,
             number_of_exposure_outliers,
             W_outcome,
             W_exposure,
             rnt_outcome,
             number_of_strata,
             linear_nonlinear_test,
             VarExp_by_Instrument_on_Exposure,
             iv_linear_nonlinear_test,
             exposure = exposure,
             outcome = outcome ) )


  out = list(strata_linear_mods = strata_linear_mods,
             strata_ivreg_mods = strata_ivreg_mods,
             summary_stats = ss_out,
             full_linear_model = lm_mod[[1]],
             null_full_gam_model = gam_mod0[[1]],
             full_gam_model = gam_mod[[1]],
             full_ivreg_model = iv_fit[[1]],
             null_full_iv_gam_model = iv_gam_mod0[[1]],
             full_iv_gam_model = iv_gam_mod[[1]],
             model_data = mod_data
             )

  return(out)

} ## end of function
