#' runs lmfit() on output of stratify_data()
#'
#' This function stratifies the working data set into either quantiles or user defined bins
#'
#' @keywords stratify
#' @param wdata a strata_data() object
#' @param outcome a single string character of the column name for the outcome or response variable
#' @param exposure a single string character of the column name for the exposure or explanatory variable
#' @param covariates a character vector that are also column names used to define variables that will be set as covariates.
#' @return returns
#' @export
#' @examples
#' stratify_lmfit()
stratify_lmfit = function( wdata, outcome, exposure, covariates){

  ## run a sapply over the wdata list object
  df_out = t( sapply(wdata, function(x){

    ### work with complete data only
    x = na.omit(x)

    ### variability check
    var_present = apply(x, 2, function(i){ length(unique(i))}) ## only column name 'strata' should have a var == 1

    if( sum(var_present == 1) == 1){
      lm_mod = lmfit( wdata = x,
                      dependent = outcome,
                      independent = exposure,
                      covariates = covariates)
    }else{
      o = rep(NA, 15); names(o) = c("n","W","Rsq","Fstat","df","dendf","beta","se","tval","P")
      lm_mod = list(summary = o)
    }

    return(lm_mod$summary)
  }) )

  ### add row names
  rownames(df_out) = paste0("strata_", 1:nrow(df_out))

  ## Return to user
  return( as.data.frame(df_out) )

}
