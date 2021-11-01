#' Identify outliers using mean and SD or median and interquartile range
#'
#' This function identifies outliers using mean and SD or median and interquartile range
#'
#' @keywords outliers mean SD median IQR
#' @param y a numeric vector from which outlier indexes will be identified.
#' @param outlier_method a single string character of "iqr" or "sd" to determine if outlier should be determined by means and sd or medians and iqr.
#' @param outlier_cutoff a single numeric value to define a cutoff value for how many iqr or sd units outlier values.
#' @return returns vector indicating the index of an outlier(s).
#' @importFrom stats sd median quantile
#' @export
#' @examples
#' id_outliers()
id_outliers = function( y , outlier_method = "iqr", outlier_cutoff = 5){
  ## verify that x is a vector.
  if(is.vector(y) == FALSE){ stop("id_outliers() parameter y is expected to be a vector of data.") }

  ## (1) using mean and SD
  meany = mean(y, na.rm = TRUE)
  sdy = sd(y, na.rm = TRUE)
  sd_cutoff = c( meany - (outlier_cutoff*sdy) , meany + (outlier_cutoff*sdy) )

  ## (2) using median and IQR
  mediany = median(y, na.rm = TRUE)
  iqry = quantile(y, probs = c(0.25, 0.75), na.rm = TRUE)
  iqry = iqry[2] - iqry[1]
  iqr_cutoff = c( mediany - (outlier_cutoff*iqry) , mediany + (outlier_cutoff*iqry) )

  ## (3) median or mean ??
  if(outlier_method == "iqr"){
    outliers = which( y < iqr_cutoff[1] | y > iqr_cutoff[2] )
  } else {
    outliers = which( y < sd_cutoff[1] | y > sd_cutoff[2] )
  }

  return(outliers)
}
