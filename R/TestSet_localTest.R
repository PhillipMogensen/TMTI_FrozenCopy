#' Test a subset of hypotheses in its closure using a user-specified local test
#'
#' @param LocalTest Function which defines a combination test.
#' @param pvals Numeric vector of p-values
#' @param subset Numeric vector; the subset to be tested
#' @param alpha Numeric; the level to test at, if stopping early. Defaults
#' to 0.05
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
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' is.sorted. Defaults to FALSE.
#' @param ... Additional arguments to be passed onto TMTI()
#'
#' @return The adjusted p-value for the test of the hypothesis that there are
#' no false hypotheses among the selected subset.
#' @export
#'
#' @examples
#' ## Simulate p-values; 10 from false hypotheses, 10 from true
#' pvals = sort(c(
#'   rbeta(10, 1, 20), # Mean value of .1
#'   runif(10)
#' ))
#' ## Test whether the highest 10 contain any false hypotheses using a Bonferroni test
#' TestSet_LocalTest(function(x) {
#'   min(c(1, length(x) * min(x)))
#' }, pvals, subset = 11:20)
TestSet_LocalTest = function(LocalTest,
                              pvals,
                              subset,
                              alpha = 0.05,
                              EarlyStop = FALSE,
                              verbose = FALSE,
                              mc.cores = 1L,
                              chunksize = 4 * mc.cores,
                              is.sorted = FALSE,
                              ...) {
  m = length(pvals)
  m2 = length(subset)
  if (is.sorted) {
    pSub = pvals[subset]
    pRest = pvals[-subset]
  } else {
    pSub = sort(pvals[subset])
    pRest = sort(pvals[-subset])
  }

  is_subset_sequence = all(seq_along(subset) == subset)

  p_first = LocalTest(pSub)
  if (p_first > alpha & EarlyStop) {
    return(p_first)
  }

  # out = list()

  # pfirst = LocalTest(pSub)
  # out[[1]] = c(
  #   "p"     = pfirst,
  #   "layer" = 0,
  #   "Accept" = (pfirst > alpha)
  # )
  #
  # if (EarlyStop & out[[1]][3]) {
  #   return(out[[1]][1])
  # }


  if (mc.cores <= 1L) {
    # stepCounter = 0
    # currentMax = 0
    #
    # out = list()
    #
    # for (i in length(pRest):1) {
    #   stepCounter = stepCounter + 1
    #   if (verbose) {
    #     cat("\rStep", stepCounter, " of ", length(pRest))
    #   }
    #
    #   p = LocalTest (
    #     if (is_subset_sequence)
    #       c(pSub, pRest[i:length(pRest)])
    #     else
    #       sort(c(pSub, pRest[i:length(pRest)]))
    #   )
    #
    #   if (p > currentMax)
    #     p = currentMax
    #
    #   if (EarlyStop & (currentMax > 0.05))
    #     break
    # }
    #
    # return(currentMax)

    out = TestSet_C(
      LocalTest = LocalTest,
      pSub = pSub,
      pRest = pRest,
      alpha = alpha,
      is_subset_sequence = is_subset_sequence,
      EarlyStop = EarlyStop,
      verbose = verbose
    )

    max(out, p_first)
  } else {
    .f = function(i) {
      if (is_subset_sequence) {
        ptilde = c(pSub, pRest[i:length(pRest)])
      } else {
        ptilde = sort(c(pSub, pRest[i:length(pRest)]))
      }
      LocalTest(ptilde)
    }
    chunks = split(rev(seq(length(pRest))), ceiling(seq(length(pRest)) / chunksize))
    results = list()
    counter = 1
    for (x in chunks) {
      if (verbose) {
        cat(sprintf("\rProcessing chunk %i of %i", counter, length(chunks)))
      }
      results_ = parallel::mclapply(
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
    max(unlist(results))
  }
}
