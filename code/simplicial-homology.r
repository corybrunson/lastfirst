library(simplextree)

# k-boundary matrix
boundary_matrix <- function(st, k) {
  stopifnot(k >= 0L, k <= st$dimension + 1L, k == as.integer(k))
  k <- as.integer(k)
  if (k == 0L)
    return(matrix(0L, nrow = 0L, ncol = st$n_simplices[[1L]]))
  if (k == st$dimension + 1L)
    return(matrix(0L, nrow = st$n_simplices[[st$dimension + 1L]], ncol = 0L))
  del_k <- straverse(
    k_simplices(st, k),
    function(tau) straverse(
      k_simplices(st, k - 1L),
      function(sigma) is_face(st, sigma, tau)
    )
  )
  class(del_k) <- "integer"
  del_k
}

# Smith normal form
transpose <- function(i, j, n) {
  stopifnot(
    n >= 1L,
    i >= 1L, i <= n,
    j >= 1L, j <= n
  )
  if (i > j) {
    ij <- c(j, i)
    i <- ij[[1L]]; j <- ij[[2L]]
    rm(ij)
  }
  c(
    if (i > 1L) seq(1L, i-1L),
    j,
    if (i < j-1L) seq(i+1L, j-1L),
    i,
    if (j < n) seq(j+1L, n)
  )
}
reduce_mod2 <- function(m, x) {
  if (x > nrow(m) || x > ncol(m)) return(m)
  # find next 1
  l <- x
  while (l <= ncol(m) && all(m[seq(x, nrow(m)), l] == 0L)) l <- l + 1L
  if (l > ncol(m)) return(m)
  k <- x - 1L + which(m[seq(x, nrow(m)), l, drop = TRUE] == 1L)[[1L]]
  # exchange rows (if necessary)
  if (k > x) m <- m[transpose(x, k, nrow(m)), , drop = FALSE]
  if (l > x) m <- m[, transpose(x, l, ncol(m)), drop = FALSE]
  # use row to clear column (if possible)
  if (x < nrow(m)) {
    i <- intersect(which(m[, x, drop = TRUE] == 1L), seq(x+1L, nrow(m)))
    m[i, ] <- sweep(m[i, , drop = FALSE], 2L, m[x, ], `+`) %% 2L
  }
  # use column to clear row (if possible)
  if (x < ncol(m)) {
    j <- intersect(which(m[x, , drop = TRUE] == 1L), seq(x+1L, ncol(m)))
    m[, j] <- sweep(m[, j, drop = FALSE], 1L, m[, x], `+`) %% 2L
  }
  # recurse
  reduce_mod2(m, x+1L)
}
smith_normal_form_mod2 <- function(m) {
  if (nrow(m) == 0L || ncol(m) == 0L) return(m)
  # ensure that x is binary
  stopifnot(all(m == m %% 2L))
  # ensure that first dimension is minimal
  tr <- nrow(m) > ncol(m)
  if (tr) m <- t(m)
  # recursive procedure
  m <- reduce_mod2(m, 1L)
  # un-transpose
  if (tr) t(m) else m
}

# cycle group ranks
hom_group_ranks_mod2 <- function(st) {
  # trivial values
  rank_data <- data.frame(
    dim = seq(0L, st$dimension),
    chains = st$n_simplices,
    cycles = c(st$n_simplices[[1L]], rep(NA_integer_, st$dimension)),
    boundaries = c(rep(NA_integer_, st$dimension), 0L)
  )
  # boundary and cycle group ranks from boundary matrices
  for (d in seq(st$dimension)) {
    m <- smith_normal_form_mod2(boundary_matrix(st, d))
    b <- as.integer(sum(diag(m)))
    rank_data$boundaries[[d]] <- b
    rank_data$cycles[[d + 1L]] <- ncol(m) - b
  }
  # calculate quotient (homology) group ranks
  transform(rank_data, homology = cycles - boundaries)
}
