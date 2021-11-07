#' summary statistics for the exposure within each strata
#'
#' This function estimates a few useful summary statistics for the exposure within each strata
#'
#' @keywords stratify summary statitics
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param exposure column name to summarize across strata
#' @return returns a data frame of summary statistics on the column of interest for each strata
#' @export
#' @examples
#' stratify_sumstats()
stratify_sumstats = function( wdata, exposure){

  df_out = t( sapply( 1:length(wdata), function(i){
    x = wdata[[i]]
    N = length(x[, exposure])
    r = range(x[, exposure], na.rm = TRUE)
    m = mean(x[, exposure], na.rm = TRUE)
    s = sd(x[, exposure], na.rm = TRUE)
    out = c(N, m, r, s); names(out) = c("N","mean","min","max", "sd")
    return(out)
  }) )
  ###
  rownames(df_out) = paste0("strata_", 1:nrow(df_out))

  ## Return to user
  return( as.data.frame(df_out) )

}
