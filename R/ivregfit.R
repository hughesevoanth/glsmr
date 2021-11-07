#' MR (ivreg) model fitting
#'
#' This function fits instrumental variable (MR, ivreg) models
#'
#' @keywords MR model
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param outcome a single string character of the column name for the outcome or dependent or response variable
#' @param exposure a single string character of the column name for the exposure or independent or explanatory variable
#' @param instrument a data frame passed to function containing necessary data for analysis
#' @param covariates a character vector that are also column names used to define variables that will be set as covariates.
#' @param weights_variable a single string character of the column name for a weights variable
#' @param rnt_outcome binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
#' @return returns a ivreg object that is a list containing (1) 'fit' a ivreg() object and (2) 'summary' a vector of summary statistics.
#' @importFrom ivreg ivreg
#' @importFrom stats formula lm residuals
#' @export
#' @examples
#' ivregfit()
ivregfit = function( wdata,
                     outcome,
                     exposure,
                     instrument,
                     covariates = NA,
                     weights_variable = NA,
                     rnt_outcome = FALSE){

  ############################################
  ## I. rank normalize the dependent|outcome ?
  ############################################
  if(rnt_outcome == TRUE){
    wdata[, outcome] = rntransform( wdata[, outcome] )
  }

  ######################
  ### II. Linear model
  ######################
  if( length(na.omit(covariates))>0 ){
    form = formula(paste0( outcome, " ~ ", paste0( covariates, collapse = " + ") ," | ", exposure, " | " , instrument ))
  } else {
    form = formula( paste0( outcome, " ~ ",  exposure, " | " , instrument ) )
  }

  #########################
  ## III. RUN the MR model
  #########################
  if( length(na.omit(weights_variable))>0) {
    iv_mod = ivreg( form, data = wdata, weights = wdata[ , weights_variable] )
  } else {
    iv_mod = ivreg(form, data = wdata )
  }

  ######################
  ### IV. Summary Stats
  ######################
  ## sample size in model
  res = residuals(iv_mod)
  iv_n = length(res); names(iv_n) = "n"

  ## normality of residuals
  W = normW(res)

  ## IV summary
  s = summary(iv_mod)
  ## model coefficients
  iv_coef = s$coefficients
  ## Report beta, se, t-value, and P-value
  iv_estimates = iv_coef[exposure,]; names(iv_estimates) = c("beta","se","tval","P")
  ## R-squared
  iv_rsq = s$r.squared; names(iv_rsq) = "rsq"
  ## Wald Test for
  iv_wald = s$waldtest[c(1,3,4,2)]
  names(iv_wald) = c("Wald_stat","Wald_df1","Wald_df2","Wald_P")
  ## Diagnostic Test
  diag_test = s$diagnostics
  ## Weak Instrument (F) Test
  iv_Ftest = diag_test[1, c(3,1,2,4)]
  names(iv_Ftest) = c("Fstat","F_df1","F_df2","F_P")
  ## Wu Hausman Endogeneity Test
  iv_Wu_Hausman = diag_test[1, c(3,2,1,4)]
  names(iv_Wu_Hausman) = c("WuH_stat","WuH_df1","WuH_df2","WuH_P")

  ######################
  ## V. Linear model data out
  ######################
  iv_out = list(fit = iv_mod,
                summary = c(iv_n, W, iv_rsq, iv_wald, iv_Ftest, iv_Wu_Hausman, iv_estimates) )

  return(iv_out)

}
