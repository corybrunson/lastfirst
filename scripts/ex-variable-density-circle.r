# packages
library(tidyverse)
library(landmark)

# source settings
source(here::here("code/settings.r"))

# initialize persistence data with parameter ranges
bumpy_persistence <- crossing(
  n = 10 ^ seq(2, 4, .5),
  p = c(0, .01, .05, .1),
  th = c(1, 5/6, 3/4, 2/3, 1/2) * pi,
  sd = c(0, 1/6, 1/4, 1/3) * pi,
  r = c(1, 2, 3, 4, 6, 10)
)

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

#' 1. Generate a sample from the bumpy circle.
#' 2. Compute towers of maxmin and lastfirst (enlarged) covers.
#' 3. Transform birth and death values to metric on S^1.
#' 4. Compute nerve filtrations of covers.
#' 5. Compute persistent homology of filtrations.
#' 6. Identify regimes where maxmin vs lastfirst better locates the 1-feature.

# loop over parameter choices
# for (i in seq_along(nrow(bumpy_persistence))) {###
i <- 583L

n <- bumpy_persistence$n[[i]]
p <- bumpy_persistence$p[[i]]
th <- bumpy_persistence$th[[i]]
sd <- bumpy_persistence$sd[[i]]
r <- bumpy_persistence$r[[i]]

#' 1. Generate a sample from the bumpy circle.

x <- bumpy_circle(n, p, th, sd, r)

#' 2. Compute towers of maxmin and lastfirst (enlarged) covers.

x_m <- landmarks_maxmin(x, pick_method = "random", cover = TRUE,
                        extend_radius = extension(mult = .1))

#' 3. Transform birth and death values to metric on S^1.



#' 4. Compute nerve filtrations of covers.



#' 5. Compute persistent homology of filtrations.



}###

#' 6. Identify regimes where maxmin vs lastfirst better locates the 1-feature.

