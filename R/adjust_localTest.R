#' Adjust all marginally significant using a Closed Testing Procedeure and a
#' user-defined local test which satisfies the quadratic shortcut given in Mogensen and Markussen (2021)
#'
#' @param localTest A function specifying a local test.
#' @param pvals vector of pvalues
#' @param alpha signicance level. Defaults to 0.05
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' is.sorted. Defaults to FALSE.
#' @param earlyStop Logical; set to TRUE to stop as soon as a hypothesis can be
#' accepted at level alpha. This speeds up the procedure, but now only provides
#' lower bounds on the p-values for the global test.
#' @param verbose Logical; set to TRUE to print progress.
#' @param mc.cores Number of cores to parallelize onto. If mc.cores > 1 the procedure
#' cannot stop early and may then be slower than with mc.cores = 1.
#' @param direction string that is equal to either "increasing"/"i" or "decreasing"/d.
#' Determines the search direction. "increasing" computes the exact adjusted p-value
#' for all those hypotheses that can be rejected (while controlling the FWER),
#' but is potentially slower than "decreasing". "decreasing" computes on the largest
#' p-value of those that can be rejected, but identifies all hypotheses that can
#' be rejected. Defaults to "increasing" and has no effect when mc.cores > 1.
#' @param ... Additional arguments
#'
#' @return a data.frame containing adjusted p-values and their respective indices.
#' @export
#'
#' @examples
#' p = sort(runif(100))  # Simulate and sort p-values
#' p[1:10] = p[1:10]**3  # Make the bottom 10 smaller, such that they correspond to false hypotheses
#' adjust_localTest (
#'     localTest = function(x) {min(c(1, length(x) * min(x)))},
#'     p, alpha = 0.05, is.sorted = TRUE
#' )
adjust_localTest = function (
    localTest,
    pvals, alpha = 0.05,
    is.sorted = FALSE,
    earlyStop = FALSE,
    verbose = FALSE,
    mc.cores = 1L,
    direction = "increasing",
    chunksize = 4 * mc.cores,
    ...
) {
  m2 = sum(pvals <= alpha)

  if (m2 <= 0)
    stop("There are no p-values that are marginally significant at level alpha")
  if (verbose)
    cat(sprintf("\rThere are %i marginally significant p-values to adjust", m2))

  if (is.sorted) {
    ord = seq_along(pvals)
  } else {
    ord <- order(pvals)
  }

  if (mc.cores == 1) {
    .f = function (i, es = earlyStop) {
      if (verbose)
        cat(sprintf("\rAdjusting p-value %i of %i.", i, m2))

      out = TestSet_localTest (
        localTest,
        pvals = pvals,
        subset = ord[i],
        alpha = alpha,
        earlyStop = es,
        verbose = verbose,
        is.sorted = is.sorted,
        ...
      )
      out
    }
    results = list()
    if (tolower(direction) == "increasing" | tolower(direction) == "i") {
      for (i in seq(m2)) {
        results[[i]] = .f(i)
        if (results[[i]] >= alpha & earlyStop) {
          message(paste0("Adjusted p-value ", i, " was above alpha, implying that the remaining are also above alpha. Exiting"))
          break
        }
      }
      return (
        data.frame (
          "p"     = unlist(results),
          "index" = ord[1:length(results)]
        )
      )
    } else if (tolower(direction) == "decreasing" | tolower(direction) == "d") {
      for (i in m2:1) {
        results[[i]] = .f(i, TRUE)
        if ((results[[i]] < alpha) & (i > 1)) {
          message(paste0("Adjusted p-value ", i, " is below alpha, implying that the remaining are also below alpha. Exiting"))
          brokeEarly = TRUE
          break
        }
      }
      if (!exists("brokeEarly") | i == 1) {
        return (
          data.frame (
            "p"     = unlist(results),
            "index" = ord[length(results):1]
          )
        )
      } else {
        nonnull = length(unlist(results))

        return (
          data.frame (
            "p"     = c (
              rep(paste0("p > ", alpha), nonnull - 1),
              results[[i]],
              rep(paste0("p <= ", results[[i]]), m2 - nonnull)
            ),
            "index" = ord[m2:1]
          )
        )
      }
    }
  } else {
    .f = function (i, es = earlyStop) {
      if (verbose)
        cat(sprintf("\rAdjusting p-value %i of %i.", i, m2))

      out = TestSet_localTest (
        localTest,
        pvals = pvals,
        subset = ord[i],
        alpha = alpha,
        earlyStop = es,
        verbose = verbose,
        is.sorted = is.sorted,
        mc.cores = mc.cores,
        chunksize = chunksize,
        ...
      )
      out
    }
    results = list()
    for (i in seq(m2)) {
      results[[i]] = .f(i)
      if (results[[i]] >= alpha & earlyStop) {
        message(paste0("Adjusted p-value ", i, " was above alpha, implying that the remaining are also above alpha. Exiting"))
        break
      }
    }
    return (
      data.frame (
        "p"     = unlist(results),
        "index" = ord[1:length(results)]
      )
    )
  }
  # else {
  #   chunks = split(seq(m2), ceiling(seq(m2) / chunksize))
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
  #       message(paste0("an adjusted p-value in chunk", counter, " was above alpha, implying that the remaining chunks are all above alpha. Exiting"))
  #       break
  #     }
  #     counter = counter + 1
  #   }
  #   return (
  #     data.frame (
  #       "p"     = unlist(results),
  #       "index" = ord[1:length(unlist(results))]
  #     )
  #   )
  # }
}
