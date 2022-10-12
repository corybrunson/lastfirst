library(tidyverse)

# source and store directories
if (dir.exists("/blue")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
  lastfirst_dir <- "~/lastfirst"
} else if (str_detect(here::here(), "jason.brunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  save_dir <- "data/cover"
  lastfirst_dir <- here::here()
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

read_rds(file.path(lastfirst_dir, str_c("data/sphere-dominance.rds"))) %>%
  ungroup() %>%
  mutate(landmarks = factor(
    landmarks,
    levels = c("maxmin", "lastfirst", "random")
  )) ->
  dominance

dominance %>%
  select(iter, distribution, duplication, landmarks, complex, relative) %>%
  ggplot(aes(y = relative)) +
  #facet_grid(complex ~ type) +
  facet_grid(complex ~ distribution + duplication) +
  geom_boxplot(aes(x = landmarks)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Relative dominance", x = "Procedure") ->
  dominance_relative_plot
# revision, with color and less annotation
dominance %>%
  select(iter, distribution, duplication, landmarks, complex, relative) %>%
  arrange(distribution, duplication, landmarks, complex) %>%
  mutate(parameters = interaction(distribution, duplication, sep = ", ")) %>%
  mutate(parameters = fct_inorder(parameters)) %>%
  ggplot(aes(y = relative)) +
  #facet_grid(complex ~ type) +
  facet_grid(complex ~ .) +
  geom_boxplot(aes(x = parameters, color = landmarks)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = proc_pal) +
  theme(legend.position = "right") +
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
# revision, with color and less annotation
dominance %>%
  select(iter, distribution, duplication, landmarks, complex, absolute) %>%
  arrange(distribution, duplication, landmarks, complex) %>%
  mutate(parameters = interaction(distribution, duplication, sep = ", ")) %>%
  mutate(parameters = fct_inorder(parameters)) %>%
  ggplot(aes(y = absolute)) +
  #facet_grid(complex ~ type) +
  facet_grid(complex ~ .) +
  geom_boxplot(aes(x = parameters, color = landmarks)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = proc_pal) +
  theme(legend.position = "right") +
  labs(y = "Absolute dominance", x = "Procedure") ->
  dominance_absolute_plot
ggsave(here::here("docs/figures/homology-sphere-absolute.pdf"),
       dominance_absolute_plot,
       width = textwidth, height = textwidth / phi)
