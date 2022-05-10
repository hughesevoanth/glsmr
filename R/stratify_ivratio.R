#' estimates the causal effect estimate or beta_iv for strata
#'
#' This function estimates the causal effect estimate or beta_iv for strata
#'
#' @keywords stratify
#' @param strata_beta_ie stratified instrument on exposure estimates from stratify_lmfit()
#' @param strata_beta_io stratified instrument on outcome estimates from stratify_lmfit()
#' @param beta_ie the instrument on exposure effect estimate for the entire data set
#' @param tsls_ss the two-stage least square summary table from the full data set
#' @return returns
#' @export
#' @examples
#' stratify_ivratio()
stratify_ivratio = function( strata_beta_ie, strata_beta_io, beta_ie, tsls_ss){
  ## Define a data frame
  df = data.frame( n =  strata_beta_ie$n )

  ### RATIO causal effect estimates using stratified beta_ie estimates
  df$beta_strata_iv = strata_beta_io$beta / strata_beta_ie$beta
  df$se_strata_iv = strata_beta_io$se / strata_beta_ie$beta

  ### RATIO causal effect SEs using the fulldata beta_ie estimate
  df$beta_iv = strata_beta_io$beta / beta_ie
  df$se_iv = strata_beta_io$se / beta_ie

  ## P-value estimates
  df$P = strata_beta_io$P

  ## ADD Exposure Summary
  x = strata_beta_ie[, c("outcome_n","outcome_mean","outcome_min","outcome_max","outcome_sd")]
  colnames(x) = gsub("outcome_","exposure_",colnames(x))
  df = cbind(df, x)

  ## ADD Outcome Summary
  x = strata_beta_io[, c("outcome_n","outcome_mean","outcome_min","outcome_max","outcome_sd")]
  df = cbind(df, x)


  ## add summary stats from full data MR to data from each strata
  temp = tsls_ss["beta_iv", c("n","beta","se","beta","se","P",
                              "exposure_n","exposure_mean","exposure_min","exposure_max","exposure_sd",
                              "outcome_n","outcome_mean","outcome_min","outcome_max","outcome_sd"
                               )]
  names(temp) = names(df)
  df = rbind(df, fulldata = temp)

  return(df)
}

