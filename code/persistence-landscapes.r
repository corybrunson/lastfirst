library(tidyverse)

phom <- as_tibble(read_rds("~/Repos/sphere-phom.rds"))
ph <- select(slice(phom, seq(16L)), dimension, birth, death)

# Algorithm 4
landscape <- function(ph, dimension, orders = Inf, tol = .Machine$double.eps) {
  
  # all orders unless otherwise specified
  stopifnot(length(orders) == 1L)
  
  # preprocess persistence data
  ph <- as.matrix(ph)
  # assume columns are dimension (optional), birth, and death
  if (ncol(ph) == 2L) {
    if (! missing(dimension))
      warning("No dimension column provided, `dimension` will be ignored.")
  } else if (ncol(ph) == 3L) {
    ph <- ph[ph[, 1L] == dimension, -1L, drop = FALSE]
  } else {
    stop("Matrix `ph` must have 2 or 3 columns: (dimension), birth, and death.")
  }
  # remove features with zero persistence
  ph <- ph[ph[, 2L] - ph[, 1L] > tol, , drop = FALSE]
  
  # Algorithm 4 (Bubenik & Dlotko, 2017)
  ph <- ph[order(ph[, 1L], -ph[, 2L]), , drop = FALSE]
  pl <- list()
  k <- 1L
  while (nrow(ph) > 0L && k <= orders) {
    plk <- matrix(NA_real_, nrow = 0L, ncol = 2L)
    b <- unname(ph[1L, 1L])
    d <- unname(ph[1L, 2L])
    ph <- ph[-1L, , drop = FALSE]
    p <- 1L
    if (b == -Inf && d == Inf) {
      plk <- rbind(plk, c(-Inf, Inf), c(Inf, Inf))
    } else if (d == Inf) {
      plk <- rbind(plk, c(-Inf, 0), c(b, 0), c(Inf, Inf))
    } else if (b == -Inf) {
      plk <- rbind(plk, c(-Inf, Inf))
    } else {
      plk <- rbind(plk, c(-Inf, 0), c(b, 0), c((b + d) / 2, (d - b) / 2))
    }
    while (plk[nrow(plk), 1L] != Inf ||
           (plk[nrow(plk), 2L] != 0 &&
            plk[nrow(plk), 2L] != Inf)) {
      if (nrow(ph) == 0L || p > nrow(ph) ||
          all(d >= ph[seq(p, nrow(ph)), 2L])) {
        plk <- rbind(plk, c(d, 0), c(Inf, 0))
      } else {
        q <- min(which(d < ph[seq(p, nrow(ph)), 2L])) + p - 1L
        b_ <- unname(ph[q, 1L])
        d_ <- unname(ph[q, 2L])
        ph <- ph[-q, , drop = FALSE]
        p <- q
        if (b_ > d) {
          plk <- rbind(plk, c(d, 0))
        }
        if (b_ >= d) {
          plk <- rbind(plk, c(b_, 0))
        } else {
          plk <- rbind(plk, c((b_ + d) / 2, (d - b_) / 2))
          ph <- rbind(ph, c(b_, d))
          ph <- rbind(
            ph[if (p == 1L) c() else seq(p - 1L), , drop = FALSE],
            ph[seq(p, nrow(ph)), , drop = FALSE][order(
              ph[seq(p, nrow(ph)), 1L], -ph[seq(p, nrow(ph)), 2L]
            ), , drop = FALSE]
          )
          p <- q + 1L
        }
        if (d_ == Inf) {
          plk <- rbind(plk, c(Inf, Inf))
        } else {
          plk <- rbind(plk, c((b_ + d_) / 2, (d_ - b_) / 2))
          b <- b_
          d <- d_
        }
      }
    }
    pl[[k]] <- plk
    k <- k + 1L
  }
  pl
}

pl_0 <- landscape(ph, dimension = 0L)
pl_1 <- landscape(ph, dimension = 1L)
pl_2 <- landscape(ph, dimension = 2L)

landscape_supp <- function(pl) {
  if (length(pl) == 0L || nrow(pl[[1L]]) <= 2L) return(double(0))
  range(pl[[1L]][seq(2L, nrow(pl[[1L]]) - 1L), 1L, drop = TRUE])
}

layer_eval <- function(p, x, tol = .Machine$double.eps) {
  # identify interval over which to compute each value
  x_cut <- cut(x, breaks = p[, 1L], include.lowest = TRUE)
  # intervals for which calculations must be done
  x_lev <- match(sort(unique(x_cut)), levels(x_cut))
  # no contribution from unsupported intervals
  x_lev <- setdiff(x_lev, c(1L, nrow(p) - 1L))
  # array of contributions
  res <- rowSums(vapply(x_lev, FUN.VALUE = x, function(l) {
    w <- which(as.integer(x_cut) == l)
    v <- rep(0, length(x))
    v[w] <- p[l, 2L] + (p[l + 1L, 2L] - p[l, 2L]) *
      (x[w] - p[l, 1L]) / (p[l + 1L, 1L] - p[l, 1L])
    v
  }))
  res[abs(res) <= tol] <- 0
  res
}

layer_fun <- function(p) {
  function(x) {
    # identify interval over which to compute each value
    x_cut <- cut(x, breaks = p[, 1L], include.lowest = TRUE)
    # intervals for which calculations must be done
    x_lev <- match(sort(unique(x_cut)), levels(x_cut))
    # no contribution from unsupported intervals
    x_lev <- setdiff(x_lev, c(1L, nrow(p) - 1L))
    # array of contributions
    rowSums(vapply(x_lev, FUN.VALUE = x, function(l) {
      w <- which(as.integer(x_cut) == l)
      v <- rep(0, length(x))
      v[w] <- p[l, 2L] + (p[l + 1L, 2L] - p[l, 2L]) *
        (x[w] - p[l, 1L]) / (p[l + 1L, 1L] - p[l, 1L])
      v
    }))
  }
}

x <- seq(0, 2, by = .01)
plot(x, layer_eval(pl_1[[1L]], x), type = "l")
f <- layer_fun(pl_1[[1]])
plot(x, f(x), type = "l")

# validation against TDA package
PH <- as.matrix(ph[! is.infinite(ph$death), , drop = FALSE])
PH_tseq <- seq(min(PH[, 2:3]), max(PH[, 2:3]), length = 500)
PH_1 <- TDA::landscape(PH, dimension = 1L, tseq = PH_tseq)
plot(x = PH_tseq, y = PH_1, type = "l")

# expand layer functionality to landscapes
landscape_eval <- function(pl, x) {
  lapply(pl, function(p) layer_eval(p, x))
}

plot(x, landscape_eval(pl_1, x)[[1L]], type = "l")
lines(x, landscape_eval(pl_1, x)[[2L]], type = "l", col = 2L)

# multiply a persistence landscape layer by a scalar
layer_scal <- function(p, a) {
  p[, 2L] <- p[, 2L] * a
  p
}
# add two persistence landscape layers
layer_sum <- function(p1, p2, tol = .Machine$double.eps) {
  x <- sort(unique(c(p1[, 1L], p2[, 1L])))
  while (any(abs(diff(x)) <= tol)) {
    x <- x[c(TRUE, abs(diff(x)) > tol)]
  }
  y <- layer_eval(p1, x) + layer_eval(p2, x)
  cbind(x, y)
}

x <- seq(0, 2, by = .01)
f <- layer_fun(layer_sum(pl_1[[1L]], pl_1[[2L]]))
plot(x, f(x), type = "l")
f <- layer_fun(layer_sum(layer_scal(pl_1[[1L]], 2), pl_1[[2L]]))
plot(x, f(x), type = "l")
f <- layer_fun(layer_sum(layer_scal(pl_1[[1L]], .5), pl_1[[2L]]))
plot(x, f(x), type = "l")

# additive identity landscape
landscape_zero <- list(matrix(c(-Inf, Inf, 0, 0), nrow = 2L, ncol = 2L))
# multiply a persistence landscape by a scalar
landscape_scal <- function(pl, a) {
  lapply(pl, function(p) layer_scal(p, a))
}
# add two persistence landscapes
landscape_sum <- function(pl1, pl2) {
  if (length(pl1) < length(pl2)) {
    pl1 <- c(pl1, replicate(length(pl2) - length(pl1), landscape_zero))
  } else if (length(pl2) < length(pl1)) {
    pl2 <- c(pl2, replicate(length(pl1) - length(pl2), landscape_zero))
  }
  mapply(layer_sum, p1 = pl1, p2 = pl2, SIMPLIFY = FALSE)
}
landscape_cumsum <- function(lst) {
  res <- landscape_zero
  for (i in seq_along(lst)) {
    res <- landscape_sum(res, lst[[i]])
  }
  res
}
landscape_mean <- function(lst) {
  res <- landscape_cumsum(lst)
  res <- landscape_scal(res, 1/length(lst))
  res
}

landscape_scal(pl_1, 2)
landscape_sum(pl_1, pl_2)
plot(x, landscape_eval(landscape_sum(pl_1, pl_2), x)[[1L]], type = "l",
     lwd = 2)
lines(x, landscape_eval(landscape_sum(pl_1, pl_2), x)[[2L]], type = "l",
      lty = 2L, col = 2L)

# distances between persistence landscapes
interval_integral <- function(x0, x1, y0, y1, p, tol = .Machine$double.eps) {
  if (abs(x0 - x1) <= tol ||
      (abs(y0) <= tol && abs(y1) <= tol))
    return(0)
  if (abs(y0 - y1) <= tol)
    return((y0 ^ p) * (x1 - x0))
  a <- (y1 - y0) / (x1 - x0)
  b <- y0 - a * x0
  abs(((a * x1 + b) ^ (p + 1) - (a * x0 + b) ^ (p + 1)) / (a * (p + 1)))
}
layer_norm <- function(plk, p) {
  ints <- vapply(
    seq(nrow(plk) - 1L),
    function(i) interval_integral(
      plk[i, 1L], plk[i + 1L, 1L], plk[i, 2L], plk[i + 1L, 2L], p
    ),
    NA_real_
  )
  sum(ints)
}
landscape_dist <- function(pl1, pl2, p = 2) {
  pl_diff <- landscape_sum(pl1, landscape_scal(pl2, -1))
  if (! is.numeric(p) || p < 1) {
    stop("Distances are only defined for norms `1 <= p <= Inf`.")
  } else if (p == Inf) {
    return(max(abs(do.call(rbind, pl_diff)[, 2L, drop = TRUE])))
  } else {
    pl_norms <- sapply(pl_diff, function(pli) layer_norm(pli, p = p))
    sum(pl_norms) ^ (1/p)
  }
}

# samples from sphere
sph_pl <- lapply(seq(100L), function(i) {
  sph <- tdaunif::sample_sphere(100L, dim = 2L)
  sph <- sph / mean(dist(sph))
  sph_ph <- ripserr::vietoris_rips(sph, dim = 2L)
  landscape(sph_ph, dimension = 1L, orders = 6L)
})
sph_pl_mean <- landscape_mean(sph_pl)
x <- landscape_supp(sph_pl_mean)
x <- seq(x[[1L]], x[[2L]], by = 0.001)
sph_pl_mean_val <- landscape_eval(sph_pl_mean[seq(6L)], x)
plot(x, sph_pl_mean_val[[1L]], type = "l", asp = 1, lwd = 2)
for (i in seq(2L, 6L)) {
  lines(x, sph_pl_mean_val[[i]], type = "l", lty = i, col = i)
}
