#' linear model fitting to estimate exposure values in the absence of the instrument effect
#'
#' This function estimates exposure values in the absence of the instrument effect by fitting a linear model of instrument + covariates on exposure and extracting the residuals.
#'
#' @keywords linear model
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param exposure a single string character of the column name for the exposure variable, which will be set as the response variable in this model.
#' @param instrument a single string character of the column name for the instrument, which will be set as the independent or explanatory variable in the model.
#' @param covariates a character vector that are also column names used to define variables that will be set as covariates.
#' @param exposure_mean_normalize binary TRUE or FALSE. TRUE statements will re-center the residulas to have a mean value of that observed for the exposure.
#' @return returns the data frame provided but with an additional column 'iv_free_exposure'.
#' @export
#' @examples
#' iv_free_exposure()
iv_free_exposure = function( wdata,
                  exposure,
                  instrument,
                  covariates = NULL,
                  exposure_mean_normalize = TRUE){

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
  res = residuals(lm_mod)

  ## accounting for NAs in model, to insure length res == nrow(wdata)
  m = match(index, names(res))
  res = res[m]; names(res) = index

  #########################
  ## IV. mean normalize the
  ##     residuals
  #########################
  if(exposure_mean_normalize == TRUE){
    exposure_mean = mean(wdata[, exposure], na.rm = TRUE)
    res = res + exposure_mean
  }

  #########################
  ## IV. add IV free exposure
  ##     to working data frame
  #########################
  wdata$iv_free_exposure = res

  #########################
  ## VI. Return to user
  #########################
  return(wdata)

}
