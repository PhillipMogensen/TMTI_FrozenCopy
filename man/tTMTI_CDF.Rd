% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tTMTI_CDF.R
\name{tTMTI_CDF}
\alias{tTMTI_CDF}
\title{Computes the analytical version of the tTMTI_infty CDF. When m>100, this should
not be used.}
\usage{
tTMTI_CDF(x, m, tau)
}
\arguments{
\item{x}{Point in which to evaluate the CDF.}

\item{m}{Number of independent tests to combine.}

\item{tau}{The truncation point of the tTMTI procedure.}
}
\value{
The probability that the test statistic is at most x assuming
independence under the global null hypothesis.
}
\description{
Computes the analytical version of the tTMTI_infty CDF. When m>100, this should
not be used.
}
\examples{
tTMTI_CDF(0.05, 100, 0.05)

}
