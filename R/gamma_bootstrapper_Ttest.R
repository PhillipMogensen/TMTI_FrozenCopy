#' Compute a list of TMTI CDFs for one- and two-sample test scenarios
#'
#' @param Y A d*m matrix of m response variables with d observations. Can
#' contain missing values in places.
#' @param X Null if one-sample, a vector with only two unique values if
#' two-sample.
#' @param n Number (or Inf) indicating what kind of minimum to consider.
#' Defaults to Inf, corresponding to the global minimum.
#' @param B Number of bootstrap replicates. Rule of thumb is to use at least
#' 10 * m
#' @param log.p Logical indicating whether to calculate p-values on log-scale.
#' Defaults to FALSE.
#' @param mc.cores Integer denoting the number of cores to use when using
#' parallelization, Defaults to 1, corresponding to single-threaded computations
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#'
#' @return A list of bootstrapped TMTI CDFs that can be used directly in the
#' TMTI_CTP function
#' @export
#'
#' @examples
#' d <- 100
#' m <- 3
#'
#' X <- sample(LETTERS[1:2], d, replace = TRUE)
#' Y <- matrix(rnorm(d * m), nrow = d, ncol = m)
#' pvalues <- apply(Y, 2, function(y) t.test(y ~ X)$p.value)
#'
#' gammaFunctions <- gamma_bootstrapper_Ttest(Y, X)  # Produces a list of CDFs
#' TMTI_CTP(pvalues, gammaList = gammaFunctions)  # Adjusted p-values using the bootstrapped CDFs

gamma_bootstrapper_Ttest <- function (
  Y,
  X = NULL,
  n = Inf,
  B = 1e3,
  log.p = FALSE,
  mc.cores = 1L,
  tau = NULL,
  K = NULL
) {
  if (!is.null(X)) {
    stopifnot (
      "X contains more than two unique values" = length(unique(X)) <= 2
    )

    .make_TMTI <- function (subset) {
      X2 <- sample(X)

      pvals <- lapply (
        subset,
        function (i) {
          stats::t.test(Y[, i] ~ X2)$p.value
        }
      )

      out <- TMTI::make_Y(pvals = unlist(pvals), tau = tau, K = K, log.p = log.p)

      out[.GetMinima(out, n = n)]
    }
  } else (
    .make_TMTI <- function (subset) {
      signs <- matrix (
        sample(c(-1, 1),
               nrow(Y) * ncol(Y),
               replace = T),
        nrow = nrow(Y),
        ncol = ncol(Y)
      )

      pvals <- lapply (
        subset,
        function (i) {
          stats::t.test(signs[, i] * Y[, i])$p.value
        }
      )
      out <- TMTI::make_Y(pvals = unlist(pvals), tau = tau, K = K, log.p = log.p)

      out[.GetMinima(out, n = n)]
    }
  )

  lapply (
    1:ncol(Y),
    function (i) {
      cat("\rComputing gamma function for level ", i, " of ", ncol(Y))
      if(i == 1)
        function (x) x
      else {
        forCDF <- unlist(parallel::mclapply (
          1:B,
          function (j) .make_TMTI(sample(1:ncol(Y), i)),
          mc.cores = mc.cores
        ))

        function (x) mean(forCDF <= x)
      }
    }
  )
}

if(F) {
  d <- 100
  m <- 3
  X <- sample(LETTERS[1:2], d, replace = T)
  Y <- matrix(rnorm(d * m), nrow = d, ncol = m)
  gamma_bootstrapper_Ttest(Y)
}
