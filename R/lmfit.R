#' linear model fitting
#'
#' This function fits observational linear models
#'
#' @keywords linear model
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param dependent a single string character of the column name for the dependent or response variable
#' @param independent a single string character of the column name for the independent or explanatory variable
#' @param covariates a character vector that are also column names used to define variables that will be set as covariates.
#' @param rnt_dependent binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
#' @return returns a lmfit object that is a list containing (1) 'fit' a lm() object and (2) 'summary' a vector of summary statistics.
#' @export
#' @examples
#' lmfit()
lmfit = function( wdata,
                  dependent = NA,
                  independent = NA,
                  covariates = NA,
                  rnt_dependent = FALSE){

  ############################################
  ## I. rank normalize the dependent|outcome ?
  ############################################
  if(rnt_dependent == TRUE){
    wdata[, dependent] = rntransform( wdata[, dependent] )
  }

  ######################
  ### II. Linear model
  ######################
  if( !is.na(covariates)[1] ){
    form = formula(paste0(dependent, " ~ ", paste0( covariates, collapse = " + ") , " + ", independent ))
  } else {
    form = formula(paste0(dependent, " ~ ", independent ))
  }

  #########################
  ## III. RUN the LINEAR MODEL
  #########################
  lm_mod = lm(form, data = wdata)

  ## sample size in model
  res = residuals(lm_mod)
  lm_n = length(res); names(lm_n) = "n"

  ## normality of residuals
  W = normW(res)

  ######################
  ### IV. lm stats out
  ######################
  ## model summary
  s = summary(lm_mod)
  ## model coefficients
  lm_coef = s$coefficients
  ## Report beta, se, t-value, and P-value
  lm_estimates = lm_coef[independent, ]; names(lm_estimates) = c("beta","se","tval","P")
  lm_rsq = s$r.squared; names(lm_rsq) = "Rsq"
  lm_F = s$fstatistic; names(lm_F) = c("Fstat","df","dendf")

  ######################
  ## V. Linear model data out
  ######################
  lm_out = list( fit = lm_mod,
                 summary = c(lm_n, W, lm_rsq, lm_F, lm_estimates)
                 )

  return(lm_out)

}
