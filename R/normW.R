#' Estimate shapiro-wilk's W statistic for normality
#'
#' This function estimates the Shapiro-Wilk's W statistic for normality
#'
#' @keywords normality Shapiro-Wilk W-statistic
#' @param y a numeric vector from which Shapiro-Wilk test of normality will be performed
#' @return an estimate of normality, the W-statistic.
#' @importFrom stats shapiro.test
#' @export
#' @examples
#' normW()
normW = function( y ){
  ## verify that x is a vector.
  if(is.vector(y) == FALSE){ stop("id_outliers() parameter y is expected to be a vector of data.") }

  ## verify that there is some variability
  if (length(unique(y)) == 1)
    stop("trait is monomorphic")
  if (length(unique(y)) == 2)
    stop("trait is binary")

  ## if there are more than 5000 observations
  ## sample down to n=5000
  if(length(y) <= 5000){
    W = shapiro.test( y )$stat; names(W) = "W"
  } else {
    W = shapiro.test( sample(y,5000) )$stat; names(W) = "W"
  }

  return(W)
}
