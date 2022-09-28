# packages
library(tidyverse)

# source settings
source(here::here("code/settings.r"))

# restore progress
bumpy_persistence <- read_rds(file = here::here("data/bumpy-betti.rds"))

# check progress
bumpy_persistence %>%
  transmute(m, done = map_lgl(betti, ~ length(.) > 1L)) %>%
  count(m, done)

#' Transform birth and death values to metric on S^1.

message("This step is not necessary in order to evaluate landmark persistence.")

#' Identify regimes where maxmin vs lastfirst better locates the 1-feature.

contig_len <- function(x) max(diff(which(c(0L, diff(x), 0L) != 1L))) - 1L

bumpy_persistence %>%
  # restrict to calculated Betti numbers
  filter(map_lgl(betti, ~ .[[1L]] != 0L)) %>%
  mutate(
    # calculate landmark range over which beta_0 = 1
    b0_lp = map(betti, ~ which(.[, 1L] == 1L)),
    # calculate landmark range over which beta_1 = 1
    b1_lp = map(betti, ~ which(.[, 2L] == 1L)),
    # calculate landmark range over which both Betti numbers are correct
    b01_lp = map(betti, ~ which(.[, 1L] == 1L & .[, 2L] == 1L))
  ) %>%
  # longest contiguous intervals (of landmark accumulation) with correct Bettis
  mutate(across(c(b0_lp, b1_lp, b01_lp), ~ map_int(., contig_len))) %>%
  print() ->
  bumpy_persistence

# check whether each parameter matters (to b_1)
bumpy_persistence %>%
  # filter(me > 0 & ae > 0) %>%
  group_by(n, m) %>%
  ggplot(aes(x = b1_lp, fill = m)) +
  facet_wrap(vars(n), nrow = 3L) +
  geom_histogram(position = "dodge")
# angle & ratio between bumps, standard deviations of bumps
bumpy_persistence %>%
  mutate(across(c(th, sd), round, digits = 2L)) %>%
  filter(b1_lp > 0L) %>%
  ggplot(aes(x = factor(n), y = b1_lp)) +
  facet_grid(rows = vars(sd), cols = vars(th, r)) +
  geom_boxplot(aes(color = m))
bumpy_persistence %>%
  mutate(across(c(th, sd), round, digits = 2L)) %>%
  # group_by(p, th, sd, r, m, me, ae) %>%
  filter(b1_lp > 0L) %>%
  ggplot(aes(x = n, y = b1_lp, group = interaction(p, th, sd, r, m, me, ae))) +
  facet_grid(rows = vars(th), cols = vars(sd, r)) +
  geom_point(aes(color = m), alpha = .25) +
  geom_line(aes(color = m), alpha = .25)

# check parameter effects on b_0 and b_1
bumpy_persistence %>%
  # filter(me > 0 & ae > 0) %>%
  group_by(n, m) %>%
  ggplot(aes(x = b01_lp, fill = m)) +
  facet_wrap(vars(n), nrow = 3L) +
  geom_histogram(position = "dodge")
# ratio between bumps and standard deviations of bumps
bumpy_persistence %>%
  mutate(across(c(th, sd), round, digits = 2L)) %>%
  ggplot(aes(x = factor(n), y = b01_lp)) +
  facet_grid(
    rows = vars(sd), cols = vars(r),
    label = function(labels) label_both(labels, sep = " = ")
  ) +
  geom_boxplot(aes(color = m))
