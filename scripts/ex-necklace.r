# packages
library(tidyverse)
library(tdaunif)
library(landmark)

# source settings
source(here::here("code/settings.r"))

# circle sampler
sample_circle_radius <- function(n, bins = 1L, sd = 0, r=1L) {
  theta <- sample_strat_segment(n, bins) * 2*pi
  res <- cbind(x = r*cos(theta), y = r*sin(theta))
  colnames(res) <- c("x","y")
  data.frame(add_noise(res, sd = sd))
}

set.seed(0L)
# sample a few points around a large circle for the "necklace" outline
bigCircle <- sample_circle_radius(20, sd=0.05, r=10)
# add points sampled from smaller circles to use as the "beads" on the necklace
dat <- rbind(
  bigCircle,
  sample_circle_radius(100, sd=0.1) %>% mutate(y=y+10),
  sample_circle_radius(100, sd=0.1) %>% mutate(x=x+10),
  sample_circle_radius(100, sd=0.1) %>% mutate(x=x-10),
  sample_circle_radius(100, sd=0.1) %>% mutate(y=y-10)
)

# plot the resulting data set
ggplot(dat, aes(x, y)) +
  coord_equal() +
  geom_point() +
  theme_minimal()

# num <- 5L works fine, but any number 6L or higher seems to cause the error
# using seed_index=1L also gives the error for landmarks_lastfirst
num <- 12L
lmk_mm <-
  landmarks_maxmin(as.matrix(dat), num = num, seed_index = "minmax", cover = TRUE)
lmk_lf <-
  landmarks_lastfirst(as.matrix(dat), num = num, seed_index = "firstlast", cover = TRUE)

dat %>%
  slice(c(lmk_mm$landmark, lmk_lf$landmark)) %>%
  mutate(procedure = c(
    rep("maxmin", nrow(lmk_mm)),
    rep("lastfirst", nrow(lmk_lf))
  )) %>%
  mutate(procedure = fct_inorder(procedure)) %>%
  mutate(name = str_c(
    "L",
    c(seq(0L, nrow(lmk_mm) - 1L), seq(0L, nrow(lmk_lf) - 1L))
  )) %>%
  print() ->
  lmks

dat %>%
  ggplot(aes(x = x, y = y)) +
  coord_equal() +
  facet_grid(~ procedure) +
  geom_point(shape = 16L, alpha = .5) +
  # geom_point(data = lmks, shape = 5L, size = 4, color = "#BA5216") +
  geom_point(data = lmks, shape = 16L, aes(color = name)) +
  geom_point(data = lmks, shape = 5L, size = 4, aes(color = name)) +
  ggrepel::geom_text_repel(data = lmks, aes(label = name, color = name)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position="none") ->
  #ggtitle("Maxmin and landmark samples from a necklace sample") ->
  lmks_plot
print(lmks_plot)

# ggsave("/Users/yara.skaf/Documents/GitHub/lastfirst/docs/figures/necklace-landmarks.pdf", lmks_plot,
ggsave(
  "docs/figures/necklace-landmarks.pdf", lmks_plot,
  width = grid::unit(textwidth, "cm"),
  height = grid::unit(textwidth / (phi * 1.1), "cm")
)
