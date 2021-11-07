#' linear model fitting to estimate the instrument on dependent effect coefficient (beta_ie) and the proportion of variance explained (eta-squared) in dependent by instrument
#'
#' This function estimates the instrument on dependent effect coefficient (beta) and the proportion of variance explained (eta-squared) in dependent by instrument
#'
#' @keywords linear model
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param dependent a single string character of the column name for the dependent or response (y) variable.
#' @param instrument a single string character of the column name for the instrument, which will be set as the independent or explanatory variable in the model.
#' @param covariates a character vector that are also column names used to define variables that will be set as covariates.
#' @return returns a lmfit object that is a list containing (1) 'fit' a lm() object and (2) 'summary' a vector of summary statistics.
#' @export
#' @examples
#' iv_estimates()
iv_estimates = function( wdata, dependent, instrument, covariates){

  wdata = na.omit(wdata)
  ####################################
  ### I. Estimate Variance Explained
  ###    Instrument on dependent
  ####################################
  u_form = formula(paste0(dependent, " ~ ", instrument ))
  u_mod = lm(u_form, data = wdata)
  s = summary(u_mod)
  VarExp = s$r.squared; names(VarExp) = c("univariate_varexp")

  ####################################
  ### II. Linear model
  ####################################
  if( length(na.omit(covariates))>0 ){
    form = formula(paste0(dependent, " ~ ", paste0( covariates, collapse = " + ") , " + ", instrument ))
  } else {
    form = formula(paste0(dependent, " ~ ", instrument ))
  }

  ####################################
  ## II. RUN LINEAR MODEL
  ####################################
  mod = lm(form, data = wdata)

  ## sample size in model
  res = residuals(mod)
  n = length(res); names(n) = "n"

  ## normality of residuals
  W = normW(res)

  ### effect estimates
  s = summary(mod)
  coef = s$coefficients[instrument, c(1,2,4) ]; names(coef) = c("beta","se","P")
  ### variance explained by model
  rsq = s$r.squared; names(rsq) = c("rsq")
  ### variance explained by instrument
  a = anova(mod)
  eta_sq = a[,2]/sum(a[,2]); names(eta_sq) = rownames(a)
  etasq = eta_sq[instrument]; names(etasq) = c("etasq")

  #########################
  ## V. Return to user
  #########################
  out = c(n, W, rsq, etasq, coef, VarExp)

  #########################
  ## V. Return to user
  #########################
  return(out)
}
