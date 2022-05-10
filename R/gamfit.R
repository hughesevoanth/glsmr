#' Generalized additive model fitting
#'
#' This function fits a GAM to the data
#'
#' @keywords GAM
#' @param wdata a data frame passed to function containing necessary data for analysis.
#' @param dependent a single string character of the column name for the dependent or response variable.
#' @param independent a single string character of the column name for the independent or explanatory variable that will be modeled as a smooth. This can be left NA, and you can instead add the exposure trait to the list of "linear_covariates", to define a null model where the exposure is modeled as a parametric term.
#' @param linear_covariates a vector of string(s) that are also column names used to define variables that will be set as parametric (linear) covariates.
#' @param smooth_covariates a vector of string(s) that are also column names used to define variables that will be set as non-linear (smooth, s()) covariates.
#' @param rnt_dependent binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
#' @param exposure_trait a single string character of the column name with the focus exposure trait of interest. If the "dependent" is set to NA, the string defined here will be used to estimate exposure summary stats (n, mean, min, max, sd) for your model.
#' @return returns a list object containing two objects, (1) 'fit' which holds the gam() model and (2) 'summary' which is a vector of summary statistics derived from the model.
#' @importFrom mgcv gam s
#' @importFrom stats na.omit formula residuals
#' @export
#' @examples
#' gamfit()
gamfit = function( wdata,
                  dependent = NA,
                  independent = NA,
                  linear_covariates = NA,
                  smooth_covariates = NA,
                  rnt_dependent = FALSE,
                  return_prediction = FALSE,
                  new_prediction_data = NA,
                  exposure_trait = NA){

  ############################################
  ### 0. Define Model Data
  ############################################
  model_variables = na.omit( c(dependent, independent, linear_covariates, smooth_covariates) )
  wdata =  wdata[, c(model_variables)]

  ############################################
  ## I. rank normalize the dependent|outcome ?
  ############################################
  if(rnt_dependent == TRUE){
    wdata[, dependent] = rntransform( wdata[, dependent] )
  }

  ######################
  ### II. GAM model
  ######################
  linear_covariates = na.omit(linear_covariates)
  smooth_covariates = na.omit(smooth_covariates)
  covariates =c(linear_covariates, smooth_covariates )

  if(is.na(independent)){
        if( length(na.omit(linear_covariates))>0 & length(na.omit(smooth_covariates))>0 ){
          form = formula( paste0(dependent, " ~ ", paste0( linear_covariates, collapse = " + ") , " + " , paste0( "s(",  smooth_covariates, ")", collapse = " + ") ) )
        } else {
          if( length(na.omit(linear_covariates))>0 & length(na.omit(smooth_covariates))==0 ){
            form = formula( paste0(dependent, " ~ ", paste0( linear_covariates, collapse = " + ") ) )
          } else {
            if( length(na.omit(linear_covariates))==0 & length(na.omit(smooth_covariates))>0 ){
              form = formula( paste0(dependent, " ~ ",  paste0( "s(",  smooth_covariates, ")", collapse = " + ") ) )
            }
          }
        }
  } else{
    if( length(covariates) == 0 ){
      form = formula(paste0(dependent, " ~ ", "s(",independent,")" ))
    } else {
      if( length(na.omit(linear_covariates))>0 & length(na.omit(smooth_covariates))>0 ){
        form = formula(paste0(dependent, " ~ ", paste0( linear_covariates, collapse = " + ") , " + " , paste0( "s(",  smooth_covariates, ")", collapse = " + ") , " + ","s(",independent,")" ))
      } else {
        if( length(na.omit(linear_covariates))>0 & length(na.omit(smooth_covariates))==0){
          form = formula(paste0(dependent, " ~ ", paste0( linear_covariates, collapse = " + ") , " + ","s(",independent,")" ))
        } else {
          if( length(na.omit(linear_covariates))==0 & length(na.omit(smooth_covariates))>0){
            form = formula(paste0(dependent, " ~ ",  paste0( "s(",  smooth_covariates, ")", collapse = " + ") , " + ","s(",independent,")" ))
          }
        }
      }
    }
  }

  #########################
  ## III. RUN the GAM
  #########################
  gam_mod = gam( form, data = wdata,  method = "REML")

  #########################
  ## IV. Summary Stats
  #########################
  ## residuals of GAM
  res = residuals(gam_mod)
  ## sample size
  mod_n = length(res); names(mod_n) = "n"
  ## normality of residuals
  W = normW(res); names(W) = "W"
  ## Summary of GAM
  s = summary(gam_mod)
  ## parametric values
  pd = s$p.table
  ## smooth values
  sd = s$s.table
  if(length(sd)>0){
    if(nrow(sd)==1){
      smooth_data = sd[1,]; names(smooth_data) = paste0(rownames(sd),"_", c("edf","df","F","P") )
    }else {
      smooth_data = c()
      for(i in 1:nrow(sd)){
        o = sd[i,]; names(o) = paste0(rownames(sd)[i], "_", c("edf","df","F","P"))
        smooth_data = c(smooth_data, o)
      }
    }
  } else {
    smooth_data = NA
  }

  ## variance explained
  rsq = s$r.sq; names(rsq) = "rsq"
  dev_exp = s$dev.expl; names(dev_exp) = "dev_exp"

  ## AIC
  aic = AIC(gam_mod); names(aic) = "AIC"
  ## log likelihood
  loglik = logLik(gam_mod); names(loglik) = "loglik"
  ## REML
  reml = summary(gam_mod)$sp.criterion[1]; names(reml) = "REML"

  ######################
  ## V. Exposure | Independent summary stats
  ######################
  if(!is.na(independent)){
    ex = gam_mod$model[, independent]
    n = length(ex)
    mean = mean(ex)
    min = min(ex)
    max = max(ex)
    sd = sd(ex)
    indepenent_stats = c(n, mean, min, max, sd)
    names(indepenent_stats) = paste0("exposure_", c("n","mean","min","max","sd") )
  } else {
    if(!is.na(exposure_trait)){
      ex = gam_mod$model[, exposure_trait]
      n = length(ex)
      mean = mean(ex)
      min = min(ex)
      max = max(ex)
      sd = sd(ex)
      indepenent_stats = c(n, mean, min, max, sd)
      names(indepenent_stats) = paste0("exposure_", c("n","mean","min","max","sd") )
    } else {
      indepenent_stats = rep(NA, 5)
      names(indepenent_stats) = paste0("exposure_", c("n","mean","min","max","sd") )
    }
  }

  ######################
  ## V. Outcome summary stats
  ######################
  ex = gam_mod$model[, dependent]
  n = length(ex)
  mean = mean(ex)
  min = min(ex)
  max = max(ex)
  sd = sd(ex)
  dependent_stats = c(n, mean, min, max, sd)
  names(dependent_stats) = paste0("outcome_", c("n","mean","min","max","sd") )



  #########################
  ## VII. GAM data out
  #########################
  if( is.na(smooth_data[1]) ){
    out = c(mod_n, W, rsq, dev_exp, reml, aic, loglik, indepenent_stats, dependent_stats)
  } else {
    out = c(mod_n, W, rsq, dev_exp, smooth_data, reml, aic, loglik, indepenent_stats, dependent_stats)
  }

  gam_out = list(fit = gam_mod,
                 summary = out
                 )
  ### Return data to user
  return(gam_out)

}
