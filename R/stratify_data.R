#' stratify data
#'
#' This function stratifies the working data set into either quantiles or user defined bins
#'
#' @keywords stratify
#' @param wdata a data frame passed to function containing necessary data for analysis
#' @param stratify_on column name to identify the variable to stratify on
#' @param strata a string or vector to define how strata should be defined.
#' @return returns
#' @export
#' @examples
#' stratify_data()
stratify_data = function( wdata,
                  stratify_on,
                  strata ){


  if(length(strata) == 1){
    q = quantile( wdata[, stratify_on], probs = seq(0, 1, 1/strata), na.rm = TRUE )
  } else {
    if(class(strata) == "numeric" & length(strata) >=3 ){
      q = strata
    } else {
      stop("please check strata parameter. Acceptable values are a single numeric value indicating the number of quantiles, or a numeric vector of at least length 3 to define strata boundries.")
    }
  }

  ## identify samples of each strata
  ## and add a strata label to the model data frame
  wdata$strata = as.factor( cut(wdata[, stratify_on], q, include.lowest=TRUE, labels=FALSE) )

  ### make a data frame for each strata
  number_of_strata = length(levels(wdata$strata))

  strata_data = lapply(1:number_of_strata, function(i){
    w = which(wdata$strata == i)
    wdata[w,]
  })

  ## add names to each strata
  names(strata_data) = paste0("strata_", 1:length(strata_data))

  ## Return to user
  return( strata_data )

}
