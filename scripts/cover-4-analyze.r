# session library
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

# evaluation data
read_rds(file.path(lastfirst_dir, "data/eval-data.rds")) %>%
  mutate(proc = fct_recode(proc, maxmin = "mm", lastfirst = "lf")) ->
  eval_data

# partition coefficient (takes value 1 on crisp partitions)
eval_data %>%
  select(careunit, n_lmks, ext_mult, proc, part_coef) %>%
  ggplot(aes(x = n_lmks, y = part_coef,
             linetype = as.factor(ext_mult), color = careunit)) +
  facet_grid(. ~ proc) +
  geom_line(size = .5) +
  scale_color_brewer(type = "qual") +
  labs(x = "Number of landmarks", y = "MPC",
       linetype = "Extension", color = "Care unit")

# area under the ROC curve (c-statistic)
eval_data %>%
  select(careunit, n_lmks, ext_mult, proc, auc) %>%
  ggplot(aes(x = n_lmks, y = auc,
             linetype = as.factor(ext_mult), color = careunit)) +
  facet_grid(. ~ proc) +
  geom_line(size = .5) +
  scale_color_brewer(type = "qual") +
  labs(x = "Number of landmarks", y = "AUROC",
       linetype = "Extension", color = "Care unit")

# partition coefficient and c-statistic
eval_data %>%
  select(careunit, n_lmks, ext_mult, proc, part_coef, auc) %>%
  pivot_longer(c(part_coef, auc), names_to = "stat", values_to = "value") %>%
  mutate(stat = fct_inorder(fct_recode(stat,
                                       MPC = "part_coef", AUROC = "auc"))) %>%
  ggplot(aes(x = n_lmks, y = value,
             linetype = as.factor(ext_mult), color = careunit)) +
  facet_grid(stat ~ proc, scales = "free_y") +
  geom_line(size = .5) +
  scale_color_brewer(type = "qual") +
  labs(x = "Number of landmarks", y = "Statistic",
       linetype = "Extension", color = "Care unit") ->
  cover_eval_plot
ggsave(here::here("docs/figures/cover-evaluate.pdf"),
       cover_eval_plot,
       width = textwidth * 1/2, height = textwidth / phi)

# Euler characteristic
eval_data %>%
  select(careunit, n_lmks, ext_mult, proc, euler) %>%
  ggplot(aes(x = n_lmks, y = euler,
             linetype = as.factor(ext_mult), color = careunit)) +
  facet_grid(. ~ proc) +
  geom_line(size = .5) +
  scale_y_sqrt() +
  scale_color_brewer(type = "qual") +
  labs(x = "Number of landmarks", y = "Euler characteristic",
       linetype = "Extension", color = "Care unit")

# simplex counts
eval_data %>%
  select(careunit, n_lmks, ext_mult, proc, simplices) %>%
  mutate(dimension = replicate(nrow(.), list(seq(0L, 2L)))) %>%
  unnest(c(dimension, simplices)) %>%
  filter(dimension > 0L) %>%
  ggplot(aes(x = n_lmks, y = simplices,
             linetype = as.factor(ext_mult), color = careunit)) +
  facet_grid(dimension ~ proc, scales = "free_y") +
  geom_line(size = .5) +
  scale_y_sqrt() +
  scale_color_brewer(type = "qual") +
  labs(x = "Number of landmarks", y = "Number of simplices",
       linetype = "Extension", color = "Care unit") ->
  cover_simp_plot
ggsave(here::here("docs/figures/cover-simplices.pdf"),
       cover_simp_plot,
       width = textwidth * 1/2, height = textwidth / phi)
