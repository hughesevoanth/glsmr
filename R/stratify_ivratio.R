#' estimates the causal effect estimate or beta_iv for strata
#'
#' This function estimates the causal effect estimate or beta_iv for strata
#'
#' @keywords stratify
#' @param strata_beta_ie stratified instrument on exposure estimates from stratify_lmfit()
#' @param strata_beta_io stratified instrument on outcome estimates from stratify_lmfit()
#' @param beta_ie the instrument on exposure effect estimate for the entire data set
#' @param strata_exp_ss stratified exposure summary statistics from stratify_sumstats()
#' @param tsls_ss the two-stage least square summary table from the full data set
#' @return returns
#' @export
#' @examples
#' stratify_ivratio()
stratify_ivratio = function( strata_beta_ie, strata_beta_io, beta_ie, strata_exp_ss, tsls_ss){
  df = data.frame( n =  strata_beta_ie$n
                   # beta_ie = strata_beta_ie$beta,
                   # se_ie = strata_beta_ie$se,
                   # P_ie = strata_beta_ie$P,
                   # beta_io = strata_beta_io$beta,
                   # se_io = strata_beta_io$se,
                   # P_io = strata_beta_io$P
                   )

  ### causal effect estimates using stratfied beta_ie estimates
  # df$beta_strata_iv = df$beta_io / df$beta_ie
  # df$se_strata_iv = df$se_io / df$beta_ie

  df$beta_strata_iv = strata_beta_io$beta / strata_beta_ie$beta
  df$se_strata_iv = strata_beta_io$se / strata_beta_ie$beta

  ### causal effect estimates using the fulldata beta_ie estimate
  df$beta_iv = strata_beta_io$beta / beta_ie
  df$se_iv = strata_beta_io$se / beta_ie

  df$P = strata_beta_io$P

  ## add summary stats for each strata
  df = cbind(df, strata_exp_ss)

  temp = tsls_ss["beta_iv", c("n","beta","se","beta","se","P","N","mean","min","max","sd")]
  names(temp) = names(iv_strata_ratio)
  df = rbind(df, fulldata = temp)

  return(df)
}

