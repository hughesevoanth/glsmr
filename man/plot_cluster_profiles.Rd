% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_cluster_profiles.R
\name{plot_cluster_profiles}
\alias{plot_cluster_profiles}
\title{Plot the outcome exposure profiles derived from the function cluster_profiles()}
\usage{
plot_cluster_profiles(
  cluster_data,
  k = 10,
  exposure_values = NA,
  exposure_name = "Exposure",
  col_count = NULL
)
}
\arguments{
\item{cluster_data}{cluster_profiles() object}

\item{k}{which k (k-means) cluster to plot?}

\item{exposure_values}{what the x-axis values should be}

\item{exposure_name}{x-axis label}

\item{col_count}{how many columns to have in plot}
}
\value{
a ggplot
}
\description{
This function makes a plot of the outcome on exposure profiles
}
\examples{
plot_cluster_profiles()
}
\keyword{GAM}
