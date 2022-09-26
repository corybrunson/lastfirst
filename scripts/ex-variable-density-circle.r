# packages
library(tidyverse)
library(landmark)

# source settings
source(here::here("code/settings.r"))

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
  unnest(ae)

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

#' 1. Generate a sample from the bumpy circle.
#' 2. Compute towers of maxmin and lastfirst (enlarged) covers.
#' 3. Transform birth and death values to metric on S^1.
#' 4. Compute nerve filtrations of covers.
#' 5. Compute persistent homology of filtrations.
#' 6. Identify regimes where maxmin vs lastfirst better locates the 1-feature.

bumpy_persistence$betti <-
  as.list(vector(mode = "integer", length = nrow(bumpy_persistence)))

# loop over parameter choices
for (i in seq(nrow(bumpy_persistence))) {###
  
  n <- bumpy_persistence$n[[i]]
  p <- bumpy_persistence$p[[i]]
  th <- bumpy_persistence$th[[i]]
  sd <- bumpy_persistence$sd[[i]]
  r <- bumpy_persistence$r[[i]]
  m <- bumpy_persistence$m[[i]]
  me <- bumpy_persistence$me[[i]]
  ae <- bumpy_persistence$ae[[i]]
  
  #' 1. Generate a sample from the bumpy circle.
  
  x <- bumpy_circle(n, p, th, sd, r)
  
  #' 2. Compute towers of maxmin and lastfirst (enlarged) covers.
  
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
  
  #' 3. Compute Betti numbers of covers.
  
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
  
}###

stop("Finish generating results.")

# restore progress
bumpy_persistence <- read_rds(file = here::here("data/bumpy-betti.rds"))

#' 4. Transform birth and death values to metric on S^1.



#' 5. Identify regimes where maxmin vs lastfirst better locates the 1-feature.

contig_len <- function(x) max(diff(which(c(0L, diff(x), 0L) != 1L))) - 1L

bumpy_persistence %>%
  # restrict to calculated Betti numbers
  filter(map_lgl(betti, ~ .[[1L]] != 0L)) %>%
  mutate(
    # calculate landmark range over which beta_0 = 1
    b0_lp = map(betti, ~ which(.[, 1L] == 1L)),
    # calculate landmark range over which beta_1 = 1
    b1_lp = map(betti, ~ which(.[, 2L] == 1L))
  ) %>%
  # longest contiguous intervals (of landmark accumulation) with correct Bettis
  mutate(across(c(b0_lp, b1_lp), ~ map_int(., contig_len))) %>%
  print() ->
  bumpy_persistence

# check whether each parameter matters
bumpy_persistence %>%
  filter(me > 0 & ae > 0) %>%
  group_by(n, m) %>%
  ggplot(aes(x = b1_lp, fill = m)) +
  geom_histogram(position = "dodge")
bumpy_persistence %>%
  mutate(across(c(th, sd), round, digits = 2L)) %>%
  filter(b1_lp > 0L) %>%
  ggplot(aes(x = n, y = b1_lp)) +
  facet_grid(rows = vars(th), cols = vars(sd, r)) +
  geom_boxplot(aes(color = m))
bumpy_persistence %>%
  mutate(across(c(th, sd), round, digits = 2L)) %>%
  # group_by(p, th, sd, r, m, me, ae) %>%
  filter(b1_lp > 0L) %>%
  ggplot(aes(x = n, y = b1_lp, group = interaction(p, th, sd, r, m, me, ae))) +
  facet_grid(rows = vars(th), cols = vars(sd, r)) +
  geom_point(aes(color = m), alpha = .25) +
  geom_line(aes(color = m), alpha = .25)
