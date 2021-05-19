library(tidyverse)
library(tdaunif)
library(landmark)

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  save_dir <- "data/cover"
  lastfirst_dir <- here::here()
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
  lastfirst_dir <- "~/lastfirst"
} else {
  stop("Cannot recognize working directory.")
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

# number of iterations of each experiment
n_iter <- 120L
# number of sampled points (after boosting)
n_pts <- 600L
# power of skewing polynomial
p_skew <- 4L
# boosting denominator
d_boost <- 6L
# number of landmarks
n_lmk <- 12L

# procedures to include in comparison
lmk_funs <- list(
  maxmin = function(x, num)
    landmarks_maxmin(x, num = num, seed_index = "minmax", cover = TRUE),
  lastfirst = function(x, num)
    landmarks_lastfirst(x, num = num, seed_index = "firstlast", cover = TRUE)
)

# functions to sample from a sphere with different biases
samp_funs <- list(
  uniform_raw = function(n) sample_sphere(n, d = 2L, sd = 0),
  uniform_boost = function(n) {
    m <- n / d_boost
    s <- sample_sphere(m, d = 2L, sd = 0)
    # iterative (cumulative) boosting
    for (i in seq(d_boost - 1L)) {
      # rejection sampling
      t <- matrix(NA_real_, nrow = 0L, ncol = 3L)
      while (nrow(t) < m) {
        u <- s[sample(nrow(s), m, replace = TRUE), , drop = FALSE]
        p <- (acos(u[, "z"]) / pi) ^ p_skew
        d <- runif(m, 0, 1)
        t <- rbind(t, u[p > d, , drop = FALSE])
      }
      t <- t[seq(m), , drop = FALSE]
      s <- rbind(s, t)
    }
    s
  },
  skewed_raw = function(n) {
    # rejection sampling
    s <- matrix(NA_real_, nrow = 0L, ncol = 3L)
    while (nrow(s) < n) {
      t <- sample_sphere(n, d = 2L, sd = 0)
      p <- (acos(t[, "z"]) / pi) ^ p_skew
      d <- runif(n, 0, 1)
      s <- rbind(s, t[p > d, , drop = FALSE])
    }
    s <- s[seq(n), , drop = FALSE]
  },
  skewed_boost = function(n) {
    m <- n / d_boost
    # rejection sampling
    s <- matrix(NA_real_, nrow = 0L, ncol = 3L)
    while (nrow(s) < m) {
      t <- sample_sphere(m, d = 2L, sd = 0)
      p <- (acos(t[, "z"]) / pi) ^ (p_skew / 2L)
      d <- runif(m, 0, 1)
      s <- rbind(s, t[p > d, , drop = FALSE])
    }
    s <- s[seq(m), , drop = FALSE]
    # iterative (cumulative) boosting
    for (i in seq(d_boost - 1L)) {
      # rejection sampling
      t <- matrix(NA_real_, nrow = 0L, ncol = 3L)
      while (nrow(t) < m) {
        u <- s[sample(nrow(s), m, replace = TRUE), , drop = FALSE]
        p <- (acos(u[, "z"]) / pi) ^ (p_skew / 2L)
        d <- runif(m, 0, 1)
        t <- rbind(t, u[p > d, , drop = FALSE])
      }
      t <- t[seq(m), , drop = FALSE]
      s <- rbind(s, t)
    }
    s
  }
)

# visual example
set.seed(900)
# uniform sample, boosted (produces most different homologies)
# also consider skewed sample without boosting (`samp_funs[[3L]]`)
sph <- samp_funs[[2L]](n_pts)
#sph <- samp_funs[[3L]](n_pts)
pairs(sph)
poles <- data.frame(x = 0, y = 0, z = c(-1, 1))

# choose landmark selection procedure (1 = maxmin, 2 = lastfirst)
i <- 2L
lmk <- lmk_funs[[i]](sph, n_lmk)
l0 <- sph[lmk$landmark[[1L]], , drop = FALSE]

# cylindrical equal-area projection
xyz2lon <- function(x, y, z) atan(x / z) + pi * (z < 0)
xyz2coslat <- function(x, y, z) sqrt(x^2 + z^2)
lon0 <- xyz2lon(l0[, "x"], l0[, "y"], l0[, "z"])
coslat0 <- xyz2coslat(l0[, "x"], l0[, "y"], l0[, "z"])

proj <- function(mat) with(as.data.frame(mat), {
  data.frame(
    x = (xyz2lon(x, y, z) - lon0),
    y = y / coslat0
  )
})

# test plot
plot(proj(sph), asp = 1, pch = 16, cex = 1, col = "#00000044")
text(proj(poles), labels = c("S", "N"), col = "blue")
points(proj(l0), cex = 2, col = "red")

# store plots to a PDF
pdfloc <- str_c(
  lastfirst_dir, "/docs/figures/sphere-", names(lmk_funs)[[i]], ".pdf"
)
pdf(pdfloc)
ran <- apply(proj(sph), 2L, range)
for (j in seq(nrow(lmk))) {
  plot.new()
  plot.window(xlim = ran[, "x"], ylim = ran[, "y"], asp = 1)
  points(proj(sph[-lmk$cover_set[[j]], , drop = FALSE]),
         pch = 16, cex = 1, col = "#00000044")
  points(proj(sph[lmk$cover_set[[j]], , drop = FALSE]),
         pch = 16, cex = 1, col = "#FF000044")
  points(proj(sph[lmk$landmark[[j]], , drop = FALSE]),
         cex = 2, col = "black")
  text(proj(poles), labels = c("S", "N"), col = "blue")
}
dev.off()
