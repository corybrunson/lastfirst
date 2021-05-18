library(tidyverse)

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

dominance <- read_rds(file.path(lastfirst_dir,
                                str_c("data/sphere-dominance.rds")))

dominance %>%
  select(iter, distribution, duplication, landmarks, complex, relative) %>%
  ggplot(aes(y = relative)) +
  #facet_grid(complex ~ type) +
  facet_grid(complex ~ distribution + duplication) +
  geom_boxplot(aes(x = landmarks)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Relative dominance", x = "Procedure") ->
  dominance_relative_plot
ggsave(here::here("docs/figures/homology-sphere-relative.pdf"),
       dominance_relative_plot,
       width = textwidth, height = textwidth / phi)

dominance %>%
  select(iter, distribution, duplication, landmarks, complex, absolute) %>%
  ggplot(aes(y = absolute)) +
  #facet_grid(complex ~ type) +
  facet_grid(complex ~ distribution + duplication) +
  geom_boxplot(aes(x = landmarks)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Absolute dominance", x = "Procedure") ->
  dominance_absolute_plot
ggsave(here::here("docs/figures/homology-sphere-absolute.pdf"),
       dominance_absolute_plot,
       width = textwidth, height = textwidth / phi)
