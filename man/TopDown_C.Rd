% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{TopDown_C}
\alias{TopDown_C}
\title{Leading NA}
\usage{
TopDown_C(LocalTest, pvals, alpha)
}
\arguments{
\item{LocalTest}{A function that returns a double in (0, 1).}

\item{pvals}{A vector of p-values.}

\item{alpha}{A double indicating the significance level}
}
\description{
Computes a confidence set for the number of false hypotheses among all hypotheses
}
