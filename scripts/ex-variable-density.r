# packages
library(tidyverse)
library(landmark)
library(patchwork)

# source settings
source(here::here("code/settings.r"))

# simulated risk model
sim <- data.frame(
  n = c(70, 20, 10),
  mean = c(0, 1, 0.5),
  sd = c(0.05, 0.05, 0.5)
)

# simulated risk data
set.seed(0L)
sim %>%
  mutate(x = pmap(list(n, mean, sd), ~ rnorm(..1, mean = ..2, sd = ..3))) %>%
  mutate(y = 0) %>%
  unnest(x) %>%
  select(x, y) ->
  dat
#dat <- cbind(
#  x = c(rnorm(70,0,0.05), rnorm(20,1,0.05), rnorm(10,0.5,0.5)),
#  y = rep(0, 100)
#)

# landmark selection
num <- 4L
lmk_mm <-
  landmarks_maxmin(as.matrix(dat), num = num, seed_index = "minmax", cover = TRUE)
lmk_lf <-
  landmarks_lastfirst(as.matrix(dat), num = num, seed_index = "firstlast", cover = TRUE)

# single landmarks data set
bind_rows(
  mutate(lmk_mm, name = str_c("L", seq(0L, num - 1L)), procedure = "maxmin"),
  mutate(lmk_lf, name = str_c("L", seq(0L, num - 1L)), procedure = "lastfirst")
) %>%
  #mutate(name = map(name, ~ parse(text = .))) %>%
  mutate(procedure = fct_inorder(procedure)) %>%
  mutate(x = map_dbl(landmark, ~ dat$x[.])) %>%
  mutate(xran = map(cover_set, ~ as.data.frame(as.list(range(dat$x[.]))))) %>%
  mutate(xran = map(xran, ~ set_names(., c("xmin", "xmax")))) %>%
  unnest(xran) ->
  lmks

# visualize distribution
sim %>%
  mutate(prop = n / sum(n)) %>%
  mutate(component = row_number()) %>%
  mutate(x = replicate(
    nrow(sim),
    seq(min(dat$x), max(dat$x), length.out = 120L),
    simplify = FALSE
  )) %>%
  mutate(y = pmap(list(x, mean, sd), ~ dnorm(..1, ..2, ..3))) %>%
  unnest(c(x, y)) %>%
  mutate(y = prop * y) %>%
  select(component, x, y) %>%
  group_by(x) %>%
  summarize(y = sum(y)) ->
  dens
ggplot(dens, aes(x, y)) +
  geom_area(alpha = .2) +
  geom_point(data = dat, y = 0, alpha = .5, shape = 16L) +
  scale_y_continuous(position = "right") +
  labs(x = NULL, y = NULL) +
  ggtitle("Maxmin and landmark covers of a variable-density risk sample") ->
  dist_plot
print(dist_plot)

# visualize covers
proc_abbr <- c(maxmin = "mm", lastfirst = "lf")
lmks %>%
  group_by(procedure) %>%
  mutate(landmark_index = factor(row_number())) %>%
  ungroup() %>%
  ggplot(aes(x = x)) +
  facet_grid(procedure ~ .) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, fill = landmark_index),
    ymin = -1, ymax = 1, alpha = .2
  ) +
  scale_fill_brewer(type = "qual", palette = 6L, guide = "none") +
  geom_point(
    data = dat, aes(x = x),
    y = 0, alpha = .5, shape = 16L
  ) +
  geom_point(aes(x = x), y = 0, size = 4, shape = 5L) +
  geom_text(aes(x = x, label = name), y = 0, vjust = -1) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()
  ) +
  labs(x = NULL, y = NULL) ->
  lmks_plot
print(lmks_plot)

# patchwork of two plots
full_plot <- dist_plot + lmks_plot + plot_layout(ncol = 1L)
print(full_plot)

ggsave(
  "docs/figures/vardens-cover.pdf", full_plot,
  width = grid::unit(textwidth, "cm"),
  height = grid::unit(textwidth / phi, "cm")
)
