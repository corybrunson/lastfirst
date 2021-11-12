#' This script samples from a highly bimodal Beta distribution and rounds the
#' results to simulate divergent risk estimates taken from a limited-precision
#' instrument (e.g. a risk score).
#'
#' The Beta distribution is bimodal because the two shape parameters a,b are
#' less than 1. As a,b approach 0, the distribution grows more bimodal, i.e. the
#' center is thinned and the ends are thickened. If a > b (respectively, a < b),
#' then more samples will be concentrated at 1 (respectively, at 0).
#'
#' The resolution p controlls the rounding (at the corresponding negative power
#' of 10).

# sample size
n <- 10000L
a <- 0.05
b <- 0.15
p <- 2

# sample from Beta distribution and round to low precision
x <- rbeta(n, a, b)
y <- round(x, digits = p)

# histogram of simulated risk scores
hist(y)

# rates of duplication of simulated risk scores (breaks at rounding resolution)
hist(y, breaks = seq(0, 1, 10^(-p)))
