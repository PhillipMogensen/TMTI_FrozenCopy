.localTest_OneTest <- function (
  localTest, k, pvals, maxStep = length(pvals), alpha = 0.05, earlyStop = F,
  verbose = F,
  ...
) {
  ord <- order(pvals)
  pvals <- sort(pvals)

  # if (is.null(gammaList))
  #   gammaList <- lapply(1:length(pvals), function(i) NULL)

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

  for (i in 2:maxStep) {
    if (verbose) cat("\rStep", round(i / m * 100, 2), "%")
    if (k < (m - i + 2))
      subP <- pvals[c(k, (m - i + 2):m)]
    else
      subP <- pvals[c((m - i + 2):m)]

    # p_TMTI      <- TMTI (
    #   subP,
    #   gamma = gammaList[[i]],
    #   B = B,
    #   log.p = log.p,
    #   tau = tau, K = K[i],
    #   ...
    # )

    p_loc <- localTest(subP)

    reject_indicator <- p_loc < alpha

    out[[i]] <- c("i" = i, "p" = p_loc, "reject" = reject_indicator)
    if (earlyStop & !reject_indicator)
      break
  }

  return (
    do.call("rbind", out)
  )
}

#' A Closed Testing Procedure for any local test using an O(n^2) shortcut.
#'
#' @param localTest A function which defines the choice of local test to use.
#' @param pvals A vector of p-values
#' @param alpha Level to perform each intersection test at. Defaults to 0.05
#' @param ... Additional arguments
#'
#' @return A data.frame containing:
#' * p_adjust: The CTP adjusted p-value, controlling the FWER strongly.
#' * Index: The original index of the unsorted p-value inputs.
#' @export localTest_CTP
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals <- c (
#'   rbeta(10, 1, 20),  ## Mean value of .05
#'   runif(10)
#' )
#' TMTI_CTP(pvals, earlyStop = TRUE)

localTest_CTP <- function (
  localTest, pvals, alpha = 0.05,
  ...
) {
  ord <- order(pvals)
  p2  <- pvals
  p   <- sort(pvals)
  m   <- length(pvals)

  Q <- matrix(0, m, m)

  Q[m, m] <- p[m]
  for (i in 1:(m - 1)) {
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
# localTest_CTP <- function (
#   localTest, pvals, alpha = 0.05,
#   earlyStop = T, verbose = T,
#   ...
# ) {
#   ord <- order(pvals)
#   pvals <- sort(pvals)
#
#   m <- length(pvals)
#
#   out <- list()
#
#   for (i in 1:m) {
#     if (verbose)
#       cat("\rStep", i)
#     CTP_i <- as.data.frame(
#       .localTest_OneTest (
#         localTest = localTest,
#         k = i,
#         pvals = pvals,
#         maxStep = m - i + 1,
#         alpha = alpha,
#         earlyStop = earlyStop,
#         ...
#       )
#     )
#     out[[i]] <- c (
#       "i" = i,
#       "p_adjust" = max(CTP_i$p),
#       "FirstAccept" = ifelse (
#         sum(CTP_i$reject == 0) >= 1,
#         min(CTP_i$i[CTP_i$reject == 0]),
#         NA
#       )
#     )
#   }
#
#   out <- as.data.frame(do.call("rbind", out))
#   out$Index <- ord
#
#   return (
#     out
#   )
# }
