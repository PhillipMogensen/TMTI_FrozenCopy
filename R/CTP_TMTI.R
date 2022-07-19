#' A Closed Testing Procedure for the TMTI using an O(n^2) shortcut
#'
#' @param pvals A vector of p-values
#' @param alpha Level to perform each intersection test at. Defaults to 0.05
#' @param B Number of bootstrap replications if gamma needs to be approximated.
#' Not used if specifying a list of functions using the gammaList argument
#' or if length(pvals) <= 100. Defaults to 1000
#' @param gammaList A list of pre-specified gamma functions. If NULL, gamma
#' functions will be approximated via bootstrap, assuming independence. Defaults
#' to NULL.
#' @param log.p Logical, indicating whether to compute Y's on log-scale.
#' Defaults to FALSE
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#' @param OnlySignificant Logical, indicating whether to only compute adjusted
#' p-values for the marginally significant p-values (TRUE) or for all observed
#' p-values (FALSE). Defaults to TRUE.
#' @param progress Logical, indicating whether or not to print progress
#' @param is.sorted Logical, indicating the p-values are pre-sorted. Defaults
#' to FALSE.
#' @param ... Additional arguments
#'
#' @return A data.frame containing:
#'
#' * p_adjust: The CTP adjusted p-value, controlling the FWER strongly.
#'
#' * Index: The original index of the unsorted p-value inputs.
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals = c(
#'   rbeta(10, 1, 20), ## Mean value of .05
#'   runif(10)
#' )
#' CTP_TMTI(pvals)
#'
CTP_TMTI = function(pvals, alpha = 0.05, B = 1e3,
                     gammaList = NULL,
                     log.p = FALSE,
                     tau = NULL, K = NULL,
                     OnlySignificant = TRUE,
                     progress = FALSE,
                     is.sorted = FALSE,
                     ...) {
  if (is.sorted) {
    ord = seq_along(pvals)
    p = pvals
  } else {
    ord = order(pvals)
    p = sort(pvals)
  }
  m = length(pvals)

  if (!is.null(K)) {
    if (length(K) < length(pvals)) {
      K = rep(K, length.out = length(pvals))
    }
  }

  Q = matrix(0, m, m)

  Q[m, m] = p[m]

  n_significant = sum(pvals <= alpha)

  if (progress) {
    count_max = if (OnlySignificant) n_significant else m
  }

  for (i in 1:(m - 1)) {
    if (progress) {
      cat(
        sprintf(
          "\rAdjusting p-value %i of %i",
          i, count_max
        )
      )
    }
    counter = m

    Q[counter, i] = p[i]

    for (j in m:(i + 1)) {
      if (OnlySignificant & i >= n_significant) {
        Q[, i] = 1
        next
      }

      counter = counter - 1

      subp = p[c(i, m:j)]
      m2 = length(subp)

      Q[counter, i] = TMTI(
        pvals = subp,
        tau = tau,
        K = K[m2],
        gamma = gammaList[[m2]],
        is.sorted = is.sorted
      )
    }
  }
  for (i in 2:m) {
    Q[1:(i - 1), i] = diag(Q)[1:(i - 1)]
  }
  adjusted_p = apply(Q, 2, max)

  data.frame(
    "p_adjusted" = adjusted_p,
    "index"      = ord
  )
}
