#' Adjust all marginally significant using a Closed Testing Procedeure and the
#' TMTI family of tests.
#'
#' @param pvals vector of pvalues
#' @param alpha signicance level. Defaults to 0.05
#' @param B Number of bootstrap replications. Only relevant if length(pvals) > 100
#' and no gammaList is supplied
#' @param gammaList A list of functions. These functions should be the CDFs of
#' the chosen TMTI test for different m.
#' @param tau Number between 0 and 1 or NULL, describing the truncation level.
#' @param K Integer between >1 and m describing the truncation index.
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' is.sorted. Defaults to FALSE.
#' @param EarlyStop Logical; set to TRUE to stop as soon as a hypothesis can be
#' accepted at level alpha. This speeds up the procedure, but now only provides
#' lower bounds on the p-values for the global test.
#' @param verbose Logical; set to TRUE to print progress.
#' @param mc.cores Number of cores to parallelize onto. If mc.cores > 1 the procedure
#' cannot stop early and may then be slower than with mc.cores = 1.
#' @param chunksize Integer indicating the size of chunks to parallelize. E.g.,
#' if setting chunksize = mc.cores, each time a parallel computation is set up,
#' each worker will perform only a single task. If mc.cores > chunksize, some
#' threads will be inactive.
#' @param ... Additional arguments
#'
#' @return a data.frame containing adjusted p-values and their respective indices.
#' @export
#'
#' @examples
#' p = sort(runif(100)) # Simulate and sort p-values
#' p[1:10] = p[1:10]**3 # Make the bottom 10 smaller, such that they correspond to false hypotheses
#' adjust_TMTI(p, alpha = 0.05, is.sorted = TRUE)
adjust_TMTI = function(pvals, alpha = 0.05, B = 1e3,
                        gammaList = NULL,
                        tau = NULL, K = NULL,
                        is.sorted = FALSE,
                        EarlyStop = FALSE,
                        verbose = FALSE,
                        mc.cores = 1L,
                        chunksize = 4 * mc.cores,
                        ...) {
  m2 = sum(pvals <= alpha)
  if (m2 <= 0) {
    stop("There are no p-values that are marginally significant at level alpha")
  }
  if (verbose) {
    cat(sprintf("\rThere are %i marginally significant p-values to adjust", m2))
  }

  if (is.sorted) {
    ord = seq_along(pvals)
  } else {
    ord = order(pvals)
  }

  .f = function(i) {
    if (verbose) {
      cat(sprintf("\rAdjusting p-value %i of %i\n", i, m2))
    }

    out = TMTI::TestSet_TMTI(
      pvals = pvals,
      subset = ord[i],
      alpha = alpha,
      EarlyStop = EarlyStop,
      tau = tau, K = K,
      gammaList = gammaList,
      verbose = verbose,
      is.sorted = is.sorted,
      mc.cores = mc.cores,
      chunksize = chunksize
    )
    out
  }

  # if (mc.cores == 1) {
  #   results = list()
  #   for (i in seq(m2)) {
  #     results[[i]] = .f(i)
  #     if (results[[i]] >= alpha & EarlyStop) {
  #       message(paste0("Adjusted p-value ", i, " was above alpha, implying that the remaining are also above alpha. Exiting"))
  #       break
  #     }
  #   }
  #   return (
  #     data.frame (
  #       "p"     = unlist(results),
  #       "index" = ord[1:length(results)]
  #     )
  #   )
  # } else {
  #   results = parallel::mclapply (
  #     seq(m2),
  #     .f,
  #     mc.cores = mc.cores
  #   )
  #   return (
  #     data.frame (
  #       "p"     = unlist(results),
  #       "index" = ord[1:length(results)]
  #     )
  #   )
  # }
}
