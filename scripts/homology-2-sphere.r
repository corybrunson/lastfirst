library(tidyverse)
library(tdaunif)
library(landmark)
library(reticulate)
use_condaenv("r-reticulate")
gd <- import("gudhi")

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
  random = function(x, num)
    sample(which(! duplicated(x)), size = num),
  maxmin = function(x, num)
    landmarks_maxmin(x, num = num, seed_index = "random"),
  lastfirst = function(x, num)
    landmarks_lastfirst(x, num = num, seed_index = "random")
)

# format GUDHI persistence data as a data frame
gudhi_persistence_df <- function(ph) {
  tbl <- data.frame(t(sapply(ph, unlist)))
  names(tbl) <- c("dimension", "birth", "death")
  tbl$dimension <- as.integer(tbl$dimension)
  tbl
}

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

# visualize four sphere samplers
sph <- data.frame()
set.seed(3530L)
for (i in seq_along(samp_funs)) {
  s <- as.data.frame(samp_funs[[i]](n_pts))
  s$sample = names(samp_funs)[[i]]
  sph <- rbind(sph, s)
}
sph %>%
  separate(sample, c("distribution", "duplication")) %>%
  mutate(distribution = fct_inorder(distribution),
         duplication = fct_inorder(duplication)) %>%
  ggplot(aes(x = x, y = z)) +
  theme_bw() +
  facet_grid(distribution ~ duplication) +
  geom_point(alpha = .25) +
  coord_equal()

# sphere
diam <- data.frame()
phom <- data.frame()
for (i in seq(n_iter)) for (j in seq_along(samp_funs)) {
  # sample points
  sph <- samp_funs[[j]](n_pts)
  #
  for (k in seq_along(lmk_funs)) {
    # generate landmarks
    lmk <- lmk_funs[[k]](sph, n_lmk)
    # append diameter data
    diam <- rbind(diam, data.frame(
      iter = i,
      sample = names(samp_funs)[[j]],
      landmarks = names(lmk_funs)[[k]],
      diameter = max(dist(sph[lmk, , drop = FALSE]))
    ))
    # simplex tree constructions
    st_const <- list(
      Rips = {
        sc <- gd$RipsComplex(points = sph[lmk, , drop = FALSE],
                             max_edge_length = 3)
        sc$create_simplex_tree(max_dimension = 3L)
      },
      alpha = {
        sc = gd$AlphaComplex(points = sph[lmk, , drop = FALSE])
        sc$create_simplex_tree(max_alpha_square = 3)
      },
      witness = {
        sc <- gd$EuclideanWitnessComplex(landmarks = sph[lmk, , drop = FALSE],
                                         witnesses = sph)
        sc$create_simplex_tree(max_alpha_square = 3, limit_dimension = 3L)
      }
    )
    #
    for (l in seq_along(st_const)) {
      # calculate persistence
      ph <- gudhi_persistence_df(st_const[[l]]$persistence())
      # append persistence data
      phom <- rbind(phom, transform(ph,
                                    iter = i,
                                    sample = names(samp_funs)[[j]],
                                    landmarks = names(lmk_funs)[[k]],
                                    complex = names(st_const)[[l]]))
    }
  }
}

phom %>%
  pivot_longer(cols = c(birth, death),
               names_to = "event", values_to = "radius") %>%
  separate(sample, c("distribution", "duplication")) %>%
  mutate(distribution = fct_inorder(distribution),
         duplication = fct_inorder(duplication),
         landmarks = fct_inorder(landmarks),
         complex = fct_inorder(complex)) %>%
  group_by(iter, distribution, duplication, landmarks, complex) %>%
  arrange(radius, dimension) %>%
  mutate(incr = ifelse(event == "birth", 1L, -1L)) %>%
  mutate(hom0 = cumsum((dimension == 0L) * incr),
         hom1 = cumsum((dimension == 1L) * incr),
         hom2 = cumsum((dimension == 2L) * incr)) %>%
  select(iter, distribution, duplication, landmarks, complex,
         radius, hom0, hom1, hom2) %>%
  group_by(iter, distribution, duplication, landmarks, complex, radius) %>%
  summarize_at(vars(starts_with("hom")), last) %>%
  print() -> betti

# de Silva & Carlsson analysis
betti %>%
  rename(r0 = "radius") %>%
  mutate(r1 = lead(r0)) %>%
  filter(r0 < r1 - 120 * .Machine$double.eps) %>%
  mutate(space = case_when(
    hom0 == 1L & hom1 == 0L & hom2 == 1L ~ "sphere",
    hom0 == 1L & hom1 == 0L & hom2 == 0L ~ "point",
    TRUE ~ NA_character_
  )) %>%
  mutate(space = fct_inorder(space)) %>%
  filter(! is.na(space)) %>%
  group_by(iter, distribution, duplication, landmarks, complex, space) %>%
  filter(space == "sphere" | r1 == Inf) %>%
  filter(r0 == max(r0)) %>%
  select(iter, distribution, duplication, landmarks, complex, r0, r1, space) %>%
  rename(`0` = r0, `1` = r1) %>%
  mutate(space = case_when(space == "sphere" ~ "R", space == "point" ~ "K")) %>%
  pivot_wider(names_from = c(space), values_from = c(`0`, `1`),
              names_glue = "{space}_{.value}") %>%
  left_join(diam %>%
              separate(sample, c("distribution", "duplication")) %>%
              mutate(distribution = fct_inorder(distribution),
                     duplication = fct_inorder(duplication),
                     landmarks = fct_inorder(landmarks)),
            by = c("iter", "distribution", "duplication", "landmarks")) %>%
  mutate(K_1 = diameter) %>% select(-diameter) %>%
  mutate(relative = (R_1 - R_0) / K_0, absolute = (R_1 - R_0) / K_1) %>%
  print() -> dominance
# greatest persistent interval
#betti %>%
#  mutate(persistence = lead(radius) - radius) %>%
#  filter(hom0 == 1L & hom1 == 0L & hom2 == 1L) %>%
#  filter(persistence == max(persistence))

write_rds(dominance,
          file.path(lastfirst_dir, str_c("data/sphere-dominance.rds")))
