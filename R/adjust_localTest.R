#' Adjust all marginally significant using a Closed Testing Procedeure and a
#' user-defined local test which satisfies the quadratic shortcut given in Mogensen and Markussen (2021)
#'
#' @param LocalTest A function specifying a local test.
#' @param pvals vector of pvalues
#' @param alpha signicance level. Defaults to 0.05
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
#' adjust_LocalTest(
#'   LocalTest = function(x) {
#'     min(c(1, length(x) * min(x)))
#'   },
#'   p, alpha = 0.05, is.sorted = TRUE
#' )
adjust_LocalTest = function(LocalTest,
                             pvals, alpha = 0.05,
                             is.sorted = FALSE,
                             EarlyStop = FALSE,
                             verbose = FALSE,
                             mc.cores = 1L,
                             chunksize = 4 * mc.cores,
                             direction = "increasing",
                             parallel.direction = "breadth",
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

  if (mc.cores == 1) {
    .f = function(i, es = EarlyStop) {
      if (verbose) {
        cat(sprintf("\rAdjusting p-value %i of %i.", i, m2))
      }

      out = TestSet_LocalTest(
        LocalTest,
        pvals = pvals,
        subset = ord[i],
        alpha = alpha,
        EarlyStop = es,
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
        if (results[[i]] >= alpha & EarlyStop) {
          message(paste0("Adjusted p-value ", i, " was above alpha, implying that the remaining are also above alpha. Exiting"))
          break
        }
      }
      return(
        data.frame(
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
        return(
          data.frame(
            "p"     = unlist(results),
            "index" = ord[length(results):1]
          )
        )
      } else {
        nonnull = length(unlist(results))

        return(
          data.frame(
            "p" = c(
              rep(paste0("p > ", alpha), nonnull - 1),
              results[[i]],
              rep(paste0("p <= ", results[[i]]), m2 - nonnull)
            ),
            "index" = ord[m2:1]
          )
        )
      }
    }
  } else if (any(tolower(parallel.direction) == c("d", "depth"))) {
    .f = function(i, es = EarlyStop) {
      if (verbose) {
        cat(sprintf("\rAdjusting p-value %i of %i.", i, m2))
      }

      out = TestSet_LocalTest(
        LocalTest,
        pvals = pvals,
        subset = ord[i],
        alpha = alpha,
        EarlyStop = es,
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
      if (results[[i]] >= alpha & EarlyStop) {
        message(paste0("Adjusted p-value ", i, " was above alpha, implying that the remaining are also above alpha. Exiting"))
        break
      }
    }
    return(
      data.frame(
        "p"     = unlist(results),
        "index" = ord[1:length(results)]
      )
    )
  } else if (any(tolower(parallel.direction) == c("b", "breadth"))) {
    chunks  = split(seq(m2), ceiling(seq(m2) / chunksize))
    results = list()
    counter = 1
    for (x in chunks) {
      results[[counter]] = unlist(parallel::mclapply (
        x,
        function (j) {
          TestSet_LocalTest (
            LocalTest,
            pvals = pvals,
            subset = ord[j],
            alpha = alpha,
            EarlyStop = EarlyStop,
            verbose = verbose,
            is.sorted = is.sorted,
            mc.cores = 1,
          )
        },
        mc.cores = mc.cores
      ))
      if ((any(unlist(results) > alpha)) & EarlyStop)
        break
    }
    return(
      data.frame(
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
