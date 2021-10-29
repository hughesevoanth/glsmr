#' Generalized additive model fitting
#'
#' This function fits a GAM to the data
#'
#' @keywords GAM
#' @param wdata a data frame passed to function containing necessary data for analysis.
#' @param dependent a single string character of the column name for the dependent or response variable.
#' @param independent a single string character of the column name for the independent or explanatory variable.
#' @param linear_covariates a vector of string(s) that are also column names used to define variables that will be set as parametric (linear) covariates.
#' @param smooth_covariates a vector of string(s) that are also column names used to define variables that will be set as non-linear (smooth, s()) covariates.
#' @param rnt_dependent binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
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
                  rnt_dependent = FALSE){

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
  n = length(res); names(n) = "n"
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

  #########################
  ## V. GAM data out
  #########################
  if( is.na(smooth_data[1]) ){
    out = c(n, W, rsq, dev_exp)
  } else {
    out = c(n, W, rsq, dev_exp, smooth_data)
  }

  gam_out = list(fit = gam_mod,
                 summary = out
                 )

  return(gam_out)

}
