#' runs ivregfit() on output of stratify_data()
#'
#' This function runs ivregfit() on a stratified data set.
#'
#' @keywords stratify
#' @param wdata a list objecft from stratify_data() made up of data frames containing necessary data for ivregfit() analysis
#' @param outcome a single string character of the column name for the outcome or dependent or response variable
#' @param exposure a single string character of the column name for the exposure or independent or explanatory variable
#' @param instrument a data frame passed to function containing necessary data for analysis
#' @param covariates a character vector that are also column names used to define variables that will be set as covariates.
#' @param weights_variable a single string character of the column name for a weights variable
#' @param rnt_outcome binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.
#' @return returns
#' @export
#' @examples
#' stratify_ivregfit()
stratify_ivregfit = function( wdata,
                              outcome,
                              exposure,
                              instrument,
                              covariates = NA,
                              weights_variable = NA,
                              rnt_outcome = FALSE){

  ## iterate over the stratified data set
  df_out = t( sapply(wdata, function(x){

    ### variability check
    var_present = apply(x, 2, function(i){ length(unique(i))})

    if( sum(var_present == 1) == 1){
      iv_mod = ivregfit( wdata = x,
                         outcome = outcome,
                         exposure = exposure,
                         instrument = instrument,
                         covariates = covariates,
                         weights_variable = NA,
                         rnt_outcome = FALSE)
      }else{
        o = rep(NA, 29); names(o) = c("n","W","rsq","Wald_stat","Wald_df1","Wald_df2","Wald_P",
                                    "Fstat","F_df1","F_df2","F_P","WuH_stat","WuH_df1","WuH_df2","WuH_P",
                                    "beta","se","tval","P",
                                    "exposure_n","exposure_mean","exposure_min","exposure_max","exposure_sd",
                                    "outcome_n","outcome_mean","outcome_min","outcome_max","outcome_sd")
        iv_mod = list(summary = o)
    }
    return(iv_mod$summary)

  }) )

  ## Return to user
  return( as.data.frame(df_out) )

}
