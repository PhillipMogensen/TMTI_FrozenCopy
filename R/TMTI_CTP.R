.TMTI_OneTest <- function (
  k, pvals, maxStep = length(pvals), alpha = 0.05, B = 1e3, earlyStop = F,
  verbose = F,
  gammaList = NULL,
  log.p = TRUE,
  tau = NULL, K = NULL,
  ...
) {
  ord <- order(pvals)
  pvals <- sort(pvals)

  if (is.null(gammaList))
    gammaList <- lapply(1:length(pvals), function(i) NULL)

  m <- length(pvals)

  out <- list()

  out[[1]] <- c("i" = 1, "p" = pvals[k], "reject" = (pvals[k] < alpha))
  if (maxStep == 1) {
    return (
      as.data.frame (
        t (
          as.matrix(out[[1]])
        )
      )
    )
  }

  if (!is.null(K)) {
    if (length(K) < length(pvals))
      K <- rep(K, length.out = length(pvals))
  }

  for (i in 2:maxStep) {
    if (verbose) cat("\rStep", round(i / m * 100, 2), "%")
    if (k < (m - i + 2))
      subP <- pvals[c(k, (m - i + 2):m)]
    else
      subP <- pvals[c((m - i + 2):m)]

    p_TMTI      <- TMTI (
      subP,
      gamma = gammaList[[i]],
      B = B,
      log.p = log.p,
      tau = tau, K = K[i],
      ...
    )

    reject_indicator <- p_TMTI < alpha

    out[[i]] <- c("i" = i, "p" = p_TMTI, "reject" = reject_indicator)
    if (earlyStop & !reject_indicator)
      break
  }

  return (
    do.call("rbind", out)
  )
}

#' A Closed Testing Procedure for the TMTI using an O(n^2) shortcut
#'
#' @param pvals A vector of p-values
#' @param alpha Level to perform each intersection test at. Defaults to 0.05
#' @param B Number of bootstrap replications when approximating gamma. Not
#' relevant if specifying a list of functions using the gammaList argument.
#' Defaults to 1000
#' @param earlyStop Logical, indicating whether to stop the search early,
#' i.e. as soon as the marginal hypothesis can not be rejected. Defaults to
#' TRUE
#' @param verbose Logical, indicating whether or not to display progress.
#' Defaults to true
#' @param gammaList A list of pre-specified gamma functions. If NULL, gamma
#' functions will be approximated via bootstrap, assuming independence. Defaults
#' to NULL.
#' @param log.p Logical, indicating whether to compute Y's on log-scale.
#' Defaults to TRUE
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#' @param ... Additional arguments
#'
#' @return A tibble containing:
#' * i: The sorted index of each p-value.
#' * p_adjust: The CTP adjusted p-value, controlling the FWER strongly.
#' * FirstAccept: The first level of the test tree at which the hypothesis could
#' not be rejected. NA if it is never rejected.
#' * Index: The original index of the unsorted p-value inputs.
#' @export TMTI_CTP
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals <- c (
#'   rbeta(10, 1, 20),  ## Mean value of .05
#'   runif(10)
#' )
#' TMTI_CTP(pvals, earlyStop = TRUE)

TMTI_CTP <- function (
  pvals, alpha = 0.05, B = 1e3,
  earlyStop = T, verbose = T,
  gammaList = NULL,
  log.p = TRUE,
  tau = NULL, K = NULL,
  ...
) {
  ord <- order(pvals)
  pvals <- sort(pvals)

  m <- length(pvals)

  out <- list()

  for (i in 1:m) {
    if (verbose)
      cat("\rStep", i)
    CTP_i <- as.data.frame(
      .TMTI_OneTest (
        k = i,
        pvals = pvals,
        maxStep = m - i + 1,
        alpha = alpha,
        B = B,
        earlyStop = earlyStop,
        gammaList = gammaList,
        log.p = log.p,
        tau = tau, K = K,
        ...
      )
    )
    out[[i]] <- c (
      "i" = i,
      "p_adjust" = max(CTP_i$p),
      "FirstAccept" = ifelse (
        sum(CTP_i$reject == 0) >= 1,
        min(CTP_i$i[CTP_i$reject == 0]),
        NA
      )
    )
  }

  out <- as.data.frame(do.call("rbind", out))
  out$Index <- ord

  return (
    out
  )
}
