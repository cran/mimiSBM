% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation.R
\name{Mat_lien_alpha}
\alias{Mat_lien_alpha}
\title{Create probality-component list for clustering per view component.}
\usage{
Mat_lien_alpha(clusters, K_barre, K)
}
\arguments{
\item{clusters}{list of link between final clustering and clustering per view component.}

\item{K_barre}{Number of clusters in the final clustering}

\item{K}{Vector of size Q, indicate the number of clusters in each component.}
}
\value{
alpha :  probality-component list for clustering per view component.
}
\description{
Create probality-component list for clustering per view component.
}
