#' Fixed effect meta analysis, heterogenity, and modulator tests of stratified estimates
#'
#' This function estimates I2, H2, and tests for residual heterogenity and mean exposure modulation in a meta analysis across exposure strata
#'
#' @keywords normality Shapiro-Wilk W-statistic
#' @param strata_table stratified summary table including beta, se, and exposure mean column
#' @param modulator_col_name a string to identify the column name that will act as the modulator. Default column names is "mean".
#' @return meta analysis summary statistics
#' @importFrom metafor rma
#' @export
#' @examples
#' meta_test()
meta_test = function( strata_table, modulator_col_name = "mean" ){
  ##############################
  ## Fixed Effect Meta Analysis
  ##############################
  ## Redefine the name of the modulator (column names)
  w = which(colnames(strata_table) == modulator_col_name)
  if(length(w) == 1){
    colnames(strata_table)[w] = "modulator"
  }

  ###
  meta_test = metafor::rma(yi = strata_table$beta,  ## Effects
                              vi = strata_table$se^2,  ## Sample variance
                              method = "FE" )             ## Fixed Effects

  ## IV heterogenity stats
  meta_het_stats = c(meta_test$I2, meta_test$H2)
  names(meta_het_stats) = c("meta_I2","meta_H2")

  ## Test of Residual Heterogeneity
  meta_het_test = c(meta_test$QE, meta_test$QEp)
  names(meta_het_test) = c("meta_het_stat","meta_het_p")

  ##############################
  ## Fixed Effect Meta with Modulator
  ## Modulator is the mean of each
  ## strat
  ##############################
  meta_mod_test <- metafor::rma(yi = strata_table$beta ~ strata_table$modulator,
                                         vi = strata_table$se^2,
                                         method="FE" )
  ## Test of Residual Heterogeneity
  meta_mod_het_test = c(meta_mod_test$QE, meta_mod_test$QEp)
  names(meta_mod_het_test) = c("meta_mod_het_stat","meta_mod_het_p")

  ## Test of Moderators (Trend test or Quadratic Test)
  ## This is the effect of the mean in each strata (the strata mean is the moderator)
  meta_strata_mean_test = c(meta_mod_test$QM, meta_mod_test$QMp)
  names(meta_strata_mean_test) = c("meta_strata_mean_stat","meta_strata_mean_p")

  ## IV meta stats
  meta_stats = c(meta_het_stats, meta_het_test, meta_mod_het_test, meta_strata_mean_test)

  ## Return data
  return( meta_stats )
}
