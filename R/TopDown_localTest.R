#' TopDown localTest algorithm for estimating a 1-alpha confidence set for the number
#' of false hypotheses among a set.
#'
#' @param localTest A function specifying a local test.
#' @param pvals A vector of p-values
#' @param subset Numeric vector specifying a subset a p-values to estimate a
#' confidence set for the number of false hypotheses for. Defaults to NULL
#' corresponding to estimating a confidence set for the number of false
#' hypotheses in the entire set.
#' @param alpha Level in [0,1] at which to generate confidence set. Defaults
#' to 0.05
#' @param verbose Logical, indicating whether or not to write out the progress.
#' Defaults to TRUE
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
#' ## Estimate the confidence set using a local Bonferroni test
#' TopDown_localTest(function(x) {min(c(1, length(x) * min(x)))}, pvals)

TopDown_localTest <- function (
  localTest,
  pvals,
  subset = NULL,
  alpha = 0.05,
  verbose = TRUE,
  mc.cores = 1L,
  chunksize = 4 * mc.cores,
  ...
) {
  ord     <- order(pvals)
  pvals   <- sort(pvals)
  m       <- length(pvals)
  t_alpha <- 0

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
      p_loc <- TestSet_localTest (
        localTest,
        pvals,
        subset2,
        alpha = alpha,
        earlyStop = TRUE,
        verbose = verbose,
        mc.cores = mc.cores,
        ...
      )
      accept <- (p_loc >= alpha)
      if (accept) {
        t_alpha <- length(subset2)
        break
      }
    }
    # if (mc.cores <= 1L) {
    #   counter <- 0
    #   top <- length(subset)
    #   for (i in seq_along(subset)) {
    #       counter <- counter + 1
    #       if (verbose) {
    #         cat(
    #           sprintf(
    #             "\rOuter step %i of %i\n", counter, length(subset)
    #           )
    #         )
    #       }
    #       subset2 <- subset[length(subset):i]
    #       p_loc <- TestSet_localTest (
    #         localTest,
    #         pvals,
    #         subset2,
    #         alpha = alpha,
    #         earlyStop = TRUE,
    #         verbose = verbose,
    #         ...
    #       )
    #       accept <- (p_loc >= alpha)
    #       if (accept) {
    #         t_alpha <- length(subset2)
    #         break
    #       }
    #     }
    #   } else {
    #   .f = function (i) {
    #     subset2 <- subset[length(subset):i]
    #     p_loc <- TestSet_localTest (
    #       localTest,
    #       pvals,
    #       subset2,
    #       alpha = alpha,
    #       earlyStop = TRUE,
    #       verbose = verbose,
    #       ...
    #     )
    #     p_loc
    #   }
    #   chunks = split(seq_along(subset), ceiling(seq_along(subset) / mc.cores))
    #   results = list()
    #   counter = 1
    #   for (x in chunks) {
    #     if (verbose)
    #       cat(sprintf("\rProcessing chunk %i of %i", counter, length(chunks)))
    #     results_ = parallel::mclapply (
    #       x,
    #       .f,
    #       mc.cores = mc.cores
    #     )
    #     results[[counter]] = unlist(results_)
    #     if (any(unlist(results) > alpha)) {
    #       break
    #     }
    #     counter = counter + 1
    #   }
    #   t_alpha = m + 1 - which(unlist(results) > alpha)[1]
    # }
  } else {
    top <- length(pvals)
    if (mc.cores <= 1) {
      for (i in 1:m) {
        if (verbose) cat("\rStep", i)
        pvals_tilde <- pvals[i:m]
        p_loc <- localTest(pvals_tilde)
        accept <- (p_loc >= alpha)
        if (accept) {
          t_alpha <- length(pvals_tilde)
          break
        }
      }
    } else {
      .f = function (i) {
        pvals_tilde <- pvals[i:m]
        p_loc <- localTest(pvals_tilde)
        p_loc
      }
      chunks = split(seq(m), ceiling(seq(m) / chunksize))
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
