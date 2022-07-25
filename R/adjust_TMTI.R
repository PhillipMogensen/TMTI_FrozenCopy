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
#' @param direction string that is equal to either "increasing"/"i" or "decreasing"/d.
#' Determines the search direction. "increasing" computes the exact adjusted p-value
#' for all those hypotheses that can be rejected (while controlling the FWER),
#' but is potentially slower than "decreasing". "decreasing" computes on the largest
#' p-value of those that can be rejected, but identifies all hypotheses that can
#' be rejected. Defaults to "increasing" and has no effect when mc.cores > 1.
#' @param parallel.direction A string that is either "breadth" or "depth"
#' (or abbreviated to "b" or "d), indicating in which direction to parallelize.
#' Breadth-first parallelization uses a more efficient C++ implementation to
#' adjust each p-value, but depth-first parallelization potentially finishes
#' faster if using early stopping (EarlyStop = TRUE) and very few hypotheses
#' can be rejected.
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
                       direction = "increasing",
                       parallel.direction = "breadth",
                       ...) {
  LocalTest = function (x) {
    TMTI::TMTI(x, tau = tau, K = K, gamma = gammaList[[length(x)]])
  }
  TMTI::adjust_LocalTest (
    LocalTest = LocalTest,
    pvals = pvals,
    alpha = alpha,
    is.sorted = is.sorted,
    EarlyStop = EarlyStop,
    verbose = verbose,
    mc.cores = mc.cores,
    chunksize = chunksize,
    direction = direction,
    parallel.direction = parallel.direction,
    ...
  )

}
