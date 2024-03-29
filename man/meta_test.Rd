% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_test.R
\name{meta_test}
\alias{meta_test}
\title{Fixed effect meta analysis, heterogenity, and modulator tests of stratified estimates}
\usage{
meta_test(strata_table, modulator_col_name = "mean")
}
\arguments{
\item{strata_table}{stratified summary table including beta, se, and exposure mean column}

\item{modulator_col_name}{a string to identify the column name that will act as the modulator. Default column names is "mean".}
}
\value{
meta analysis summary statistics
}
\description{
This function estimates I2, H2, and tests for residual heterogenity and mean exposure modulation in a meta analysis across exposure strata
}
\examples{
meta_test()
}
\keyword{Shapiro-Wilk}
\keyword{W-statistic}
\keyword{normality}
