#' TopDown TMTI algorithm for estimating a 1-alpha confidence set for the number
#' of false hypotheses among a set.
#'
#' @param pvals A vector of p-values
#' @param subset Numeric vector specifying a subset a p-values to estimate a
#' confidence set for the number of false hypotheses for. Defaults to NULL
#' corresponding to estimating a confidence set for the number of false
#' hypotheses in the entire set.
#' @param alpha Level in [0,1] at which to generate confidence set. Defaults
#' to 0.05
#' @param gammaList List of pre-specified gamma functions. If NULL, the
#' functions will be approximated by bootstrap assuming independence. Defaults
#' to NULL
#' @param verbose Logical, indicating whether or not to write out the progress.
#' Defaults to TRUE
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' is.sorted. Defaults to FALSE.
#' @param mc.cores Integer specifying the number of cores to parallelize onto.
#' @param ... Additional parameters
#'
#' @return A lower 1-alpha bound for the number of false hypotheses among the
#' set of supplied p-values
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals <- c (
#'   rbeta(10, 1, 20),  ## Mean value of .05
#'   runif(10)
#' )
#' TopDown_TMTI(pvals)

TopDown_TMTI <- function (
  pvals,
  subset = NULL,
  alpha = 0.05,
  gammaList = NULL,
  verbose = TRUE,
  tau = NULL, K = NULL,
  is.sorted = FALSE,
  mc.cores = 1L,
  ...
) {
  if (is.sorted) {
    ord = seq_along(pvals)
    p   = pvals
  } else {
    ord <- order(pvals)
    p   <- sort(pvals)
  }
  m       <- length(pvals)
  t_alpha <- 0

  if(!is.null(K))
    if(length(K) == 1)
      K <- rep(K, length(pvals))

  if (!is.null(subset) & length(subset) < length(pvals)) {
    counter <- 0
    top <- length(subset)
    for (i in seq_along(subset)) {
      counter <- counter + 1
      if (verbose) {
        cat(
          sprintf(
            "\rOuter step %i of %i\n", counter, length(subset)
          )
        )
      }
      # subset2 <- subset[length(subset):i]
      subset2 <- subset[i:length(subset)]
      p_TMTI <- TestSet_TMTI (
        pvals,
        subset2,
        alpha = alpha,
        tau = tau,
        K = K,
        earlyStop = TRUE,
        gammalist = gammaList,
        verbose = verbose,
        is.sorted = is.sorted,
        ...
      )
      accept <- (p_TMTI >= alpha)
      if (accept) {
        t_alpha <- length(subset2)
        break
      }
    }
  } else {
    if (mc.cores <= 1) {
      if (is.null(gammaList))
        gammaList <- lapply(1:length(pvals), function(i) NULL)
      top <- length(pvals)
      for (i in 1:m) {
        if (verbose) cat("\rStep", i)
        pvals_tilde <- pvals[i:m]
        p_TMTI <- TMTI(pvals_tilde, gamma = gammaList[[i]], tau = tau, K = K[i], is.sorted = is.sorted, ...)
        accept <- (p_TMTI >= alpha)
        if (accept) {
          t_alpha <- length(pvals_tilde)
          break
        }
      }
    } else {
      top <- length(pvals)
      .f = function (i) {
        pvals_tilde <- pvals[i:m]
        TMTI(pvals_tilde, gamma = gammaList[[i]], tau = tau, K = K[i], is.sorted = is.sorted, ...)
      }
      chunks = split(seq(m), ceiling(seq(m) / mc.cores))
      results = list()
      counter = 1
      for (x in chunks) {
        if (verbose)
          cat(sprintf("\rProcessing chunk %i of %i", counter, length(chunks)))
        results_ = parallel::mclapply (
          x,
          .f,
          mc.cores = mc.cores
        )
        results[[counter]] = unlist(results_)
        if (any(unlist(results) > alpha)) {
          break
        }
        counter = counter + 1
      }
      t_alpha = m + 1 - which(unlist(results) > alpha)[1]
    }

  }

  if (verbose)
    cat (
      paste0 (
        "Confidence set for the number of false hypotheses is {",
        top - t_alpha,
        ",..., ",
        top,
        "}\n"
      )
    )
  return(top - t_alpha)
}
