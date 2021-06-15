#' Computes the analytical version of the rtMTI_\infty CDF. When m>100, this should
#' not be used.
#'
#' @param x Point in which to evaluate the CDF
#' @param m Number of independent tests to combine
#'
#' @return
#' @export
#'
#' @examples
#' rtTMTI_CDF(0.05, 100, 10)
rtTMTI_CDF <- function (x, m, K) { ## This is the explicit form of gamma
  if (K >= m) {
    return(TMTI_CDF(x, m))
  }

  P <- function(x, a) { ## Constructs the necessary polynomials
    m  <- length(a)

    sum(1 / factorial(m:1) * x^(m:1) * a)
  }

  xs <- numeric(m)
  xs[1:K] <- qbeta(x, 1:K, m + 1 - 1:K)
  xs[(K + 1):m] <- xs[K]

  PP <- list()
  PP[1] <- xs[1]
  for (i in 2:m) {
    PP[[i]] <- P(xs[i], c(1, -do.call("c", PP[1:(i - 1)])))
  }

  1 - factorial(m) * (
    P(1, c(1, -do.call("c", PP[1:(m - 1)]))) -
      PP[[m]]
  )
}


if (F) {
  m <- 50
  K <- 1

  gg <- TMTI::gamma_bootstrapper(m, B = 1e5, K = K, tau = NULL, log.p = F, mc.cores = 8)
  curve(Vectorize(function (x) gg(x))(x))
  curve(Vectorize(function (x) rtTMTI_CDF(x, m, K))(x), add = T, col = 2)

}
