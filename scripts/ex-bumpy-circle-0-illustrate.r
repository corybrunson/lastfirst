# packages
library(tidyverse)
library(landmark)
library(ggtda)
library(ggforce)
library(patchwork)

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

# number of points
n <- 36L
# number of landmarks
l <- 12L
# multiplicative extension
me <- 2

set.seed(736009L)
x <- bumpy_circle(n = n, p = .05, th = 3/4 * pi, sd = 1/4 * pi, r = 3)
plot(x, pch = 16L, asp = 1)

x_mm <- landmarks_maxmin(
  x, pick_method = "random", seed_index = "random", cover = TRUE, tower = TRUE,
  num = l, extend_radius = extension(mult = me, add = .1)
)
x_lf <- landmarks_lastfirst(
  x, pick_method = "random", seed_index = "random", cover = TRUE, tower = TRUE,
  num = l, extend_cardinality = extension(mult = me, add = n/12)
)

as.data.frame(x) %>%
  mutate(landmark = row_number()) ->
  x_df
x_mm %>%
  select(-tower_max) %>%
  left_join(x_df, by = "landmark") %>%
  mutate(num = row_number()) ->
  x_mm_df
# recover radii of neighborhoods
x %>%
  proxy::dist() %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(p1 = row_number()) %>%
  pivot_longer(cols = -p1, names_to = "p2", values_to = "distance") %>%
  mutate(p2 = as.integer(p2)) ->
  x_dist
x_lf %>%
  select(-tower_max) %>%
  left_join(x_df, by = "landmark") %>%
  unnest(cover_set) %>%
  left_join(x_dist, by = c(landmark = "p1", cover_set = "p2")) %>%
  group_by(landmark, x, y) %>%
  summarize(radius = max(distance)) %>%
  ungroup() %>%
  mutate(num = row_number()) ->
  x_lf_df

# colors
face_col <- "lightgrey"

ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = .5)) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  geom_disk(
    data = x_mm_df, aes(x = x, y = y),
    radius = (1 + me) * x_mm$radius[[l]],
    fill = face_col, alpha = .15, color = NA
  ) +
  geom_point(data = as.data.frame(x), aes(x = x, y = y)) +
  geom_point(
    data = x_mm_df, aes(x = x, y = y),
    color = proc_pal[[1L]], size = 4, shape = 1L
  ) +
  geom_text(
    data = x_mm_df, aes(x = x, y = y, label = num),
    hjust = "outward", vjust = "outward"
  ) +
  labs(x = NULL, y = NULL) +
  ggtitle("Maxmin cover") ->
  plot_maxmin

ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = .5)) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  geom_circle(
    data = x_lf_df, aes(x0 = x, y0 = y, r = radius),
    fill = face_col, alpha = .15, color = NA
  ) +
  geom_point(data = as.data.frame(x), aes(x = x, y = y)) +
  geom_point(
    data = x_lf_df, aes(x = x, y = y),
    color = proc_pal[[2L]], size = 4, shape = 1L
  ) +
  geom_text(
    data = x_lf_df, aes(x = x, y = y, label = num),
    hjust = "outward", vjust = "outward"
  ) +
  labs(x = NULL, y = NULL) +
  ggtitle("Lastfirst cover") ->
  plot_lastfirst

plot_maxmin + plot_lastfirst

ggsave(
  here::here("docs/figures/bumpy-covers.pdf"),
  plot_maxmin + plot_lastfirst, width = textwidth, height = textwidth / 2
)
