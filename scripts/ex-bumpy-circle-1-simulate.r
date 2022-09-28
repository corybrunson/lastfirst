# packages
library(tidyverse)
library(landmark)

# source settings
source(here::here("code/settings.r"))

# sampler
bumpy_circle <- function(
    n,
    p = .1, # share sampled uniformly; 0 < p < .5
    th = pi, # angle between bumps; 0 < th < pi
    sd = pi/4, # standard deviations of bumps; 0 < sd < pi/2
    r = 1 # ratio between bump densities; 1 < r < 10
) {
  # partition sample between bumps
  ns <- n * c(p, (1 - p) * c(1, 1/r) / (1 + 1/r))
  # sample from each bump on the topological circle
  x <- c(
    runif(ns[[1L]], 0, 2*pi),
    rnorm(ns[[2L]], 0, sd),
    rnorm(ns[[3L]], th, sd)
  ) %% (2*pi)
  # map to the standard embedding
  cbind(x = cos(x), y = sin(x))
}

# calculate Betti numbers from maximal simplices
maxs_to_betti <- function(maxs) {
  st <- simplextree::simplex_tree(maxs)
  stpy <- interplex::as_py_gudhi_simplextree(st)
  if (st$dimension == 1) {
    stpy$insert(st$n_simplices[[1L]] + seq(3L))
  }
  stpy$compute_persistence()
  b <- stpy$betti_numbers()
  if (st$dimension == 1) {
    b[[1L]] <- b[[1L]] - 1L
  }
  b
}

# initialize persistence data with parameter ranges
bumpy_persistence <- crossing(
  # data parameters
  n = c(12L, 36L, 60L),
  p = c(0, .01, .05, .1),
  th = c(1, 3/4, 1/2) * pi,
  sd = c(1/6, 1/4, 1/3) * pi,
  r = c(1, 3, 10),
  # landmark parameters
  m = fct_inorder(c("maxmin", "lastfirst")),
  me = c(0, 1, 2)
) %>%
  rowwise() %>%
  mutate(ae = case_when(
    m == "maxmin" ~ list(c(0, .1, .2)),
    m == "lastfirst" ~ list(c(0, n / 12, n / 6))
  )) %>%
  unnest(ae) %>%
  mutate(betti = as.list(vector(
    mode = "integer",
    length = nrow(bumpy_persistence)
  )))

# loop over parameter choices
for (i in seq(nrow(bumpy_persistence))) {
  
  n <- bumpy_persistence$n[[i]]
  p <- bumpy_persistence$p[[i]]
  th <- bumpy_persistence$th[[i]]
  sd <- bumpy_persistence$sd[[i]]
  r <- bumpy_persistence$r[[i]]
  m <- bumpy_persistence$m[[i]]
  me <- bumpy_persistence$me[[i]]
  ae <- bumpy_persistence$ae[[i]]
  
  #' Generate a sample from the bumpy circle.
  
  x <- bumpy_circle(n, p, th, sd, r)
  
  #' Compute towers of maxmin and lastfirst (enlarged) covers.
  
  x_l <- switch(
    as.character(m),
    maxmin = landmarks_maxmin(
      x, pick_method = "random", seed_index = "random", tower = TRUE,
      extend_radius = extension(mult = me, add = ae)
    ),
    lastfirst = landmarks_lastfirst(
      x, pick_method = "random", seed_index = "random", tower = TRUE,
      extend_cardinality = extension(mult = me, add = ae)
    )
  )
  
  #' Compute Betti numbers of covers.
  
  x_b <- matrix(NA_integer_, nrow = nrow(x_l), ncol = 2L)
  for (j in seq(nrow(x_l))) {
    if (all(vapply(x_l$tower_max[[j]], length, 0L) <= 1L)) {
      x_b[j, ] <- c(0L, 0L)
    } else {
      bets <- maxs_to_betti(x_l$tower_max[[j]])
      x_b[j, ] <- bets
    }
  }
  
  bumpy_persistence$betti[[i]] <- x_b
  
  # save progress
  write_rds(bumpy_persistence, file = here::here("data/bumpy-betti.rds"))
  
}
