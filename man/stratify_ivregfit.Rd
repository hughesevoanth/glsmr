% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stratify_ivregfit.R
\name{stratify_ivregfit}
\alias{stratify_ivregfit}
\title{runs ivregfit() on output of stratify_data()}
\usage{
stratify_ivregfit(
  wdata,
  outcome,
  exposure,
  instrument,
  covariates = NA,
  weights_variable = NA,
  rnt_outcome = FALSE
)
}
\arguments{
\item{wdata}{a list objecft from stratify_data() made up of data frames containing necessary data for ivregfit() analysis}

\item{outcome}{a single string character of the column name for the outcome or dependent or response variable}

\item{exposure}{a single string character of the column name for the exposure or independent or explanatory variable}

\item{instrument}{a data frame passed to function containing necessary data for analysis}

\item{covariates}{a character vector that are also column names used to define variables that will be set as covariates.}

\item{weights_variable}{a single string character of the column name for a weights variable}

\item{rnt_outcome}{binary TRUE or FALSE if the dependent or response variable should be rank normal transformed.}
}
\value{
returns
}
\description{
This function runs ivregfit() on a stratified data set.
}
\examples{
stratify_ivregfit()
}
\keyword{stratify}
