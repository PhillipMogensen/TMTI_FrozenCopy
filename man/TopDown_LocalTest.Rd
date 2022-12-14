% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TopDown_LocalTest.R
\name{TopDown_LocalTest}
\alias{TopDown_LocalTest}
\alias{TopDown_localTest}
\title{TopDown LocalTest algorithm for estimating a 1-alpha confidence set for the number
of false hypotheses among a set.}
\usage{
TopDown_LocalTest(
  LocalTest,
  pvals,
  subset = NULL,
  alpha = 0.05,
  verbose = FALSE,
  mc.cores = 1L,
  chunksize = 4 * mc.cores,
  direction = "binary",
  ...
)

TopDown_localTest(
  localTest,
  pvals,
  subset = NULL,
  alpha = 0.05,
  verbose = TRUE,
  mc.cores = 1L,
  chunksize = 4 * mc.cores,
  ...
)
}
\arguments{
\item{LocalTest}{A function specifying a local test.}

\item{pvals}{A vector of p-values.}

\item{subset}{Numeric vector specifying a subset a p-values to estimate a
confidence set for the number of false hypotheses for. Defaults to NULL
corresponding to estimating a confidence set for the number of false
hypotheses in the entire set.}

\item{alpha}{Level in [0,1] at which to generate confidence set. Defaults
to 0.05.}

\item{verbose}{Logical, indicating whether or not to write out the progress.
Defaults to TRUE.}

\item{mc.cores}{Integer specifying the number of cores to parallelize onto.}

\item{chunksize}{Integer indicating the size of chunks to parallelize. E.g.,
if setting chunksize = mc.cores, each time a parallel computation is set up,
each worker will perform only a single task. If mc.cores > chunksize, some
threads will be inactive.}

\item{direction}{A string indicating whether to perform a binary search ('binary'/'b')
or decreasing ('decreasing'/'d') search. Defaults to 'binary', which has better
computational complexity.}

\item{...}{Additional parameters.}

\item{localTest}{A function specifying a local test (deprecated).}
}
\value{
A 1-alpha bound lower for the number of false hypotheses among the
specified subset of the supplied p-values
}
\description{
TopDown LocalTest algorithm for estimating a 1-alpha confidence set for the number
of false hypotheses among a set.
}
\examples{
## Simulate some p-values
## The first 10 are from false hypotheses, the next 10 are from true
pvals = c(
  rbeta(10, 1, 20), ## Mean value of .05
  runif(10)
)
## Estimate the confidence set using a local Bonferroni test
TopDown_LocalTest(function(x) {
  min(c(1, length(x) * min(x)))
}, pvals)

}
