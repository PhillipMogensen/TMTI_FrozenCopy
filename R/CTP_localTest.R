#' A Closed Testing Procedure for any local test satisfying the conditions of Mogensen and Markussen (2021) using an O(n^2) shortcut.
#'
#' @param localTest A function which defines the choice of local test to use.
#' @param pvals A vector of p-values
#' @param alpha Level to perform each intersection test at. Defaults to 0.05
#' @param OnlySignificant Logical, indicating whether to only compute adjusted
#' p-values for the marginally significant p-values (TRUE) or for all observed
#' p-values (FALSE). Defaults to TRUE.
#' @param ... Additional arguments
#'
#' @return A data.frame containing:
#' * p_adjust: The CTP adjusted p-value, controlling the FWER strongly.
#' * Index: The original index of the unsorted p-value inputs.
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals <- c (
#'   rbeta(10, 1, 20),  ## Mean value of .05
#'   runif(10)
#' )
#' ## Perform the CTP using a local Bonferroni test
#' CTP_localTest(function(x) {min(c(length(x) * min(x), 1))}, pvals)

CTP_localTest <- function (
  localTest, pvals, alpha = 0.05,
  OnlySignificant = TRUE,
  ...
) {
  ord <- order(pvals)
  p2  <- pvals
  p   <- sort(pvals)
  m   <- length(pvals)

  Q <- matrix(0, m, m)

  Q[m, m] <- p[m]

  n_significant = sum(pvals <= alpha)

  for (i in 1:(m - 1)) {
    if (OnlySignificant & i >= n_significant) {
      Q[, i] = 1
      next
    }

    counter <- m

    Q[counter, i] <- p[i]

    for (j in m:(i + 1)) {
      counter <- counter - 1

      subp <- p[c(i, m:j)]
      m2   <- length(subp)

      Q[counter, i] <- localTest(subp)
    }
  }
  for (i in 1:(m - 1)) {
    Q[i, (i + 1):m] <- diag(Q)[i]
  }
  adjusted_p <- apply(Q, 2, max)

  data.frame (
    "p_adjusted" = adjusted_p,
    "index"      = ord
  )
}
