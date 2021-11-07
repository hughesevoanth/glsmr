#' linear model fitting to estimate the instrument predicted exposure, sometimes referred to as d.hat in an MR framework.
#'
#' This function estimates instrument predicted exposure by fitting a linear model of instrument + covariates on exposure and extracting the fitted values.
#'
#' @keywords linear model
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param exposure a single string character of the column name for the exposure variable, which will be set as the response variable in this model.
#' @param instrument a single string character of the column name for the instrument, which will be set as the independent or explanatory variable in the model.
#' @param covariates a character vector that are also column names used to define variables that will be set as covariates.
#' @return returns a lmfit object that is a list containing (1) 'fit' a lm() object and (2) 'summary' a vector of summary statistics.
#' @export
#' @examples
#' iv_predicted_exposure()
iv_predicted_exposure = function( wdata,
                  exposure,
                  instrument,
                  covariates = NULL){

  ######################
  ### I. Linear model
  ######################
  if( length(na.omit(covariates))>0 ){
    form = formula(paste0(exposure, " ~ ", paste0( covariates, collapse = " + ") , " + ", instrument ))
  } else {
    form = formula(paste0(exposure, " ~ ", instrument ))
  }

  #########################
  ## II. RUN the LINEAR MODEL
  #########################
  index = rownames(wdata)
  lm_mod = lm(form, data = wdata)

  #########################
  ## III. Extract the residuals
  #########################
  ## IV free exposures a.k.a the residuals
  fv = fitted.values(lm_mod)

  ## accounting for NAs in model, to insure length res == nrow(wdata)
  m = match(index, names(fv))
  fv = fv[m]; names(fv) = index

  #########################
  ## IV. add IV free exposure
  ##     to working data frame
  #########################
  wdata$iv_predicted_exposure = fv

  #########################
  ## VI. Return to user
  #########################
  return(wdata)

}
