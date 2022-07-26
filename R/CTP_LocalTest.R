#' A Closed Testing Procedure for any local test satisfying the conditions of Mogensen and Markussen (2021) using an O(n^2) shortcut.
#'
#' @name CTP_LocalTest
#' @aliases localTest_CTP
#' @param LocalTest A function which defines the choice of local test to use.
#' @param pvals A vector of p-values
#' @param alpha Level to perform each intersection test at. Defaults to 0.05
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' is.sorted. Defaults to FALSE.
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
#' pvals = c(
#'   rbeta(10, 1, 20), ## Mean value of .05
#'   runif(10)
#' )
#' ## Perform the CTP using a local Bonferroni test
#' CTP_LocalTest(function(x) {
#'   min(c(length(x) * min(x), 1))
#' }, pvals)
#'

CTP_LocalTest = function(LocalTest, pvals, alpha = 0.05, is.sorted = FALSE, ...) {
  if (is.sorted) {
    ord = 1:length(pvals)
  } else {
    ord = order(pvals)
    pvals = sort(pvals)
  }
  f = function (x, y) {
    TMTI::TestSet_C (
      LocalTest = LocalTest,
      pSub = x,
      pRest = y,
      alpha = 0.05,
      is_subset_sequence = TRUE,
      EarlyStop = FALSE,
      verbose = FALSE
    )
  }

  p_adjusted = FullCTP_C (
    LocalTest,
    f,
    pvals
  )
  data.frame (
    "p_adjusted" = p_adjusted,
    "index"      = ord
  )
}

#'
#' @rdname CTP_LocalTest
#' @param localTest A function specifying a local test.
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals = c(
#'   rbeta(10, 1, 20), ## Mean value of .05
#'   runif(10)
#' )
#' ## Perform the CTP using a local Bonferroni test
#' CTP_LocalTest(function(x) {
#'   min(c(length(x) * min(x), 1))
#' }, pvals)
#'

localTest_CTP = function(localTest, pvals, alpha = 0.05, is.sorted = FALSE, ...) {
  if (is.sorted) {
    ord = 1:length(pvals)
  } else {
    ord = order(pvals)
    pvals = sort(pvals)
  }
  f = function (x, y) {
    TMTI::TestSet_C (
      localTest = localTest,
      pSub = x,
      pRest = y,
      alpha = 0.05,
      is_subset_sequence = TRUE,
      EarlyStop = FALSE,
      verbose = FALSE
    )
  }
  .Deprecated(new = "CTP_LocalTest")

  p_adjusted = FullCTP_C (
    localTest,
    f,
    pvals
  )
  data.frame (
    "p_adjusted" = p_adjusted,
    "index"      = ord
  )
}
