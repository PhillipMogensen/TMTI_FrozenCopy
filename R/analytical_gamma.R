library(tidyverse)

##### The P/gg functions WORK -- but only for m < 135 #####
P <- function(x, a) { ## Constructs the necessary polynomials
  m  <- length(a)
  #m2 <- m - 1
  # 1/factorial(m) * x^(m) * a[1] - sum(1 / factorial(1:m2) * x^(1:m2) * a[2:m])

  sum(1 / factorial(m:1) * x^(m:1) * a)
  # sum(1 / gmp::factorialZ(m:1) * x^(m:1) * a)
}
logP <- function (x, a) {
  m  <- length(a)
  #m2 <- m - 1
  # 1/factorial(m) * x^(m) * a[1] - sum(1 / factorial(1:m2) * x^(1:m2) * a[2:m])

  sum (
    exp (
      -lfactorial(m:1) + (m:1)*log(x)
    ) * a
  )
  # exp(lfactorial(m) + m * log(x) + a[1]) - sum (
  #   exp(lfactorial((m-1):1) + ((m-1):1) * log(x) + a[2:m])
  # )
}
gg <- Vectorize(function (x, m) { ## This is the explicit form of gamma
  xs <- qbeta(x, 1:m, m + 1 - 1:m)

  PP <- list()
  PP[1] <- xs[1]
  for (i in 2:m) {
    PP[[i]] <- P(xs[i], c(1, -do.call("c", PP[1:(i - 1)])))
  }

  1 - factorial(m) * (
    P(1, c(1, -do.call("c", PP[1:(m - 1)]))) -
    PP[[m]]
  )
}, vectorize.args = "x")

gg_restricted <- Vectorize(function (x, m, maxterm = m) { ## This is the explicit form of gamma
  xs <- qbeta(x, 1:m, m + 1 - 1:m)

  PP <- list()
  PP[1] <- xs[1]
  for (i in 2:(m - 1)) {
    PP[[i]] <- logP(xs[i], c(1, -do.call("c", PP[1:(i - 1)])))
  }

  if (maxterm >= m) {
    1 - (
      1 - xs[m]^m - sum (
        # gmp::factorialZ(m) / gmp::factorialZ((m-1):1) * do.call("c", PP) * (1 - xs[m]^((m-1):1))
        # factorial(m) / factorial((m-1):1) * do.call("c", PP) * (1 - xs[m]^((m-1):1))
        exp(lfactorial(m) - lfactorial((m-1):1)) * do.call("c", PP) * (1 - xs[m]^((m-1):1))
      )
    )
  } else {
    mm <- m - maxterm
    1 - (
      1 - xs[m]^m - sum (
        # factorial(m) / factorial((m-1):mm) * do.call("c", PP)[1:mm] * (1 - xs[m]^((m-1):mm))
        exp(lfactorial(m) - lfactorial((m-1):1)) * do.call("c", PP)[1:mm] * (1 - xs[m]^((m-1):mm))
      )
    )
  }
}, vectorize.args = "x")

P2 <- function (x, a) {
  m <- length(a)

  vec <- sapply (
    m:1,
    function (i) {
      # prod(i:1)
      sum(log(i:1))
    }
  )
  # sum (
  #   # 1 / factorial(m:1) * x^(m:1) * a
  #   exp (
  #     -vec + (m:1) * log(x) + log(-a)
  #   )
  # )
  #exp(-vec[1] + m * log(x) + log(a[1]))  - sum(-vec[2:m] + (m:2) * log(x) + log(a[2:m]))

  # exp(-vec[1] + (m) * log(x) + log(a[1])) - sum(exp(-vec[2:m] + ((m-1):1) * log(x) + log(a[2:m])))

  # exp(-vec[1] + (m) * log(x))*a[1] - sum(exp(-vec[2:m] + ((m-1):1) * log(x))*a[2:m])

  sum(exp(-vec + (m:1) * log(x))*a)
}
gg2 <- Vectorize(function (x, m) { ## This is the explicit form of gamma
  xs <- qbeta(x, 1:m, m + 1 - 1:m)

  PP <- list()
  PP[1] <- xs[1]
  for (i in 2:m) {
    # print(PP[[i-1]])
    PP[[i]] <- P2(xs[i], c(1, do.call("c", PP[1:(i - 1)])))
  }

  # 1 - factorial(m) * (
  #   P(1, c(1, -do.call("c", PP[1:(m - 1)]))) -
  #     PP[[m]]
  # )
  1 - exp (
    sum(log(m:1)) +
    log(P2(1, c(1, do.call("c", PP[1:(m - 1)]))) - PP[[m]])
  )
}, vectorize.args = "x")

P3 <- function (x, a, m2) {
  m <- length(a)

  vec <- sapply (
    m:1,
    function (i) {
      # prod(i:1)
      sum(log(m2:1)) -
      sum(log(i:1))
    }
  )
  sum (
    # 1 / factorial(m:1) * x^(m:1) * a
    exp (
      vec + (m:1) * log(x) + log(a)
    )
  )
}
gg3 <- Vectorize(function (x, m) { ## This is the explicit form of gamma
  xs <- qbeta(x, 1:m, m + 1 - 1:m)

  PP <- list()
  PP[1] <- xs[1]
  for (i in 2:m) {
    PP[[i]] <- P(xs[i], c(1, -do.call("c", PP[1:(i - 1)])))
  }

  1 - exp (
    sum(log(m:1)) +
      log(P(1, c(1, -do.call("c", PP[1:(m - 1)]))) - PP[[m]])
  )
}, vectorize.args = "x")

gg4 <- Vectorize(function (x, m) {
  xs <- qbeta(x, 1:m, m + 1 - 1:m)

  PP <- list()
  PP[1] <- xs[1]
  for (i in 2:(m - 1)) {
    PP[[i]] <- P2(xs[i], c(1, -do.call("c", PP[1:(i - 1)])))
  }

  vec  <- 2:m
  vec2 <- (m + 1 - vec)
  vec3 <- cumsum(log(m:2))

  x^m + sum (
    (1 - x^vec2) * do.call("c", PP) * exp(vec3)
  )
})

gg5 <- Vectorize(function (x, m) {
  xs <- qbeta(x, 1:m, m + 1 - 1:m)

  PP <- list()
  PP[1] <- xs[1]
  for (i in 2:(m - 1)) {
    PP[[i]] <- P2(xs[i], c(1, -do.call("c", PP[1:(i - 1)])))
  }

  return(do.call("c", PP))
  return (
    factorial(m) / factorial((m-1):1) * do.call("c", PP) * (1 - xs[m]^((m-1):1))
  )
})



### This is just for checking
m <- 1000
forCDF <- lapply (
  1:1e4,
  function (i) min(pbeta(sort(runif(m)), 1:m, m + 1 - 1:m))
) %>% do.call("c", .)
gammaSim <- Vectorize(function(x) mean(forCDF <= x))

curve(gammaSim(x))
curve(Vectorize(function(x) gg(x, m))(x), add = T, col = 2, lty = 2)
curve(Vectorize(function(x) gg2(x, m))(x), add = T, col = 3, lty = 2)



curve(Vectorize(function(x) gg (x, 1e2))(x))
curve(Vectorize(function(x) gg2(x, 1e2))(x), add = T, col = 2, lty = 2)
