% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obs_modeling.R
\name{obs_modeling}
\alias{obs_modeling}
\title{linear and non-linear (GAM) observational modeling}
\usage{
obs_modeling(
  wdata,
  outcome = NA,
  exposure = NA,
  linear_covariates = NA,
  smooth_covariates = NA,
  rnt_outcome = FALSE,
  weights_variable = NA,
  outlier_method = "iqr",
  outlier_cutoff = 5,
  messages = FALSE
)
}
\arguments{
\item{wdata}{a data frame passed to function containing necessary data for analysis}

\item{outcome}{a single string character of the column name for the outcome or dependent or response variable}

\item{exposure}{a single string character of the column name for the exposure or independent or explanatory variable}

\item{linear_covariates}{a vector of string(s) that are also column names used to define variables that will be set as parametric (linear) covariates.}

\item{smooth_covariates}{a vector of string(s) that are also column names used to define variables that will be set as non-linear (smooth, s()) covariates.}

\item{rnt_outcome}{binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.}

\item{weights_variable}{a single string character of the column name for a weights variable}

\item{outlier_method}{a single string character of "iqr" or "sd" to determine if outlier should be determined by means and sd or medians and iqr.}

\item{outlier_cutoff}{a single numeric value to define a cutoff value for how many iqr or sd units outlier values}

\item{messages}{should a progress message be printed to screen - binary TRUE or FALSE}
}
\value{
returns a obs_modeling vector of summary statistics
}
\description{
This performs traditional linear and GAM or non-linear modeling
}
\examples{
obs_modeling()
}
\keyword{GAM}
\keyword{lm}
\keyword{modeling}
\keyword{observational}
