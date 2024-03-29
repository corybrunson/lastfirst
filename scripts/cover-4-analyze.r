# session library
library(tidyverse)

# source and store directories
if (dir.exists("/blue")) {
  # HiPerGator
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
  lastfirst_dir <- "~/lastfirst"
} else if (str_detect(here::here(), "jason.brunson")) {
  # laptop
  save_dir <- "data/cover"
  lastfirst_dir <- here::here()
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

# MIMIC-III

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
  scale_color_manual(values = unit_pal) +
  labs(x = "Number of landmarks", y = "Crispness (MPC)",
       linetype = "Extension", color = "Care unit")

# area under the ROC curve (c-statistic)
eval_data %>%
  select(careunit, n_lmks, ext_mult, proc, auc) %>%
  ggplot(aes(x = n_lmks, y = auc,
             linetype = as.factor(ext_mult), color = careunit)) +
  facet_grid(. ~ proc) +
  geom_line(size = .5) +
  scale_color_manual(values = unit_pal) +
  labs(x = "Number of landmarks", y = "Accuracy (AUROC)",
       linetype = "Extension", color = "Care unit")

# partition coefficient and c-statistic
eval_data %>%
  select(careunit, n_lmks, ext_mult, proc, part_coef, auc) %>%
  pivot_longer(c(part_coef, auc), names_to = "stat", values_to = "value") %>%
  mutate(stat = fct_inorder(fct_recode(
    stat,
    `Crispness (MPC)` = "part_coef",
    `Accuracy (AUROC)` = "auc"
  ))) %>%
  ggplot(aes(x = n_lmks, y = value,
             linetype = as.factor(ext_mult), color = careunit)) +
  facet_grid(stat ~ proc, scales = "free_y") +
  geom_line(size = .5) +
  scale_color_manual(values = unit_pal) +
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
  scale_color_manual(values = unit_pal) +
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
  scale_color_manual(values = unit_pal) +
  labs(x = "Number of landmarks", y = "Number of simplices",
       linetype = "Extension", color = "Care unit") ->
  cover_simp_plot
ggsave(here::here("docs/figures/cover-simplices.pdf"),
       cover_simp_plot,
       width = textwidth * 1/2, height = textwidth / phi)

# Mexican covid-19 data

# evaluation data
read_rds(file.path(lastfirst_dir, "data/eval-data-mx.rds")) %>%
  mutate(proc = fct_recode(proc, maxmin = "mm", lastfirst = "lf")) ->
  eval_data_mx

# Euler characteristic
eval_data_mx %>%
  select(n_lmks, ext_mult, proc, euler) %>%
  ggplot(aes(x = n_lmks, y = euler, linetype = as.factor(ext_mult))) +
  facet_grid(. ~ proc) +
  geom_line(size = .5) +
  scale_y_sqrt() +
  labs(x = "Number of landmarks", y = "Euler characteristic",
       linetype = "Extension")

# simplex counts
eval_data_mx %>%
  select(n_lmks, ext_mult, proc, simplices) %>%
  mutate(dimension = replicate(nrow(.), list(seq(0L, 2L)))) %>%
  unnest(c(dimension, simplices)) %>%
  filter(dimension > 0L) %>%
  ggplot(aes(x = n_lmks, y = simplices, linetype = as.factor(ext_mult))) +
  facet_grid(dimension ~ proc, scales = "free_y") +
  geom_line(size = .5) +
  scale_y_sqrt() +
  labs(x = "Number of landmarks", y = "Number of simplices",
       linetype = "Extension") ->
  cover_simp_plot_mx
ggsave(here::here("docs/figures/cover-simplices-mx.pdf"),
       cover_simp_plot_mx,
       width = textwidth * 1/2, height = textwidth / phi)

# partition coefficient (takes value 1 on crisp partitions)
eval_data_mx %>%
  select(n_lmks, ext_mult, proc, part_coef) %>%
  ggplot(aes(x = n_lmks, y = part_coef, linetype = as.factor(ext_mult))) +
  facet_grid(. ~ proc) +
  geom_line(size = .5) +
  labs(x = "Number of landmarks", y = "Crispness (MPC)",
       linetype = "Extension")

# area under the ROC curve (c-statistic)
eval_data_mx %>%
  pivot_longer(cols = c(auc_intubado, auc_def),
               names_to = "outcome", values_to = "auc") %>%
  mutate(outcome = str_remove(outcome, "^auc_")) %>%
  mutate(outcome = fct_inorder(outcome)) %>%
  mutate(outcome = fct_recode(outcome,
                              intubation = "intubado", death = "def")) %>%
  select(n_lmks, ext_mult, proc, outcome, auc) %>%
  ggplot(aes(x = n_lmks, y = auc,
             linetype = as.factor(ext_mult), color = outcome)) +
  facet_grid(. ~ proc) +
  geom_line(size = .5) +
  labs(x = "Number of landmarks", y = "Accuracy (AUROC)",
       linetype = "Extension", color = "Outcome")

# partition coefficient and c-statistic
eval_data_mx %>%
  select(n_lmks, ext_mult, proc, part_coef, auc_intubado, auc_def) %>%
  pivot_longer(c(part_coef, auc_intubado, auc_def),
               names_to = "stat", values_to = "value") %>%
  mutate(stat = fct_inorder(fct_recode(
    stat,
    `Crispness (MPC)` = "part_coef",
    `AUROC, Intubation` = "auc_intubado",
    `AUROC, Mortality` = "auc_def"
  ))) %>%
  ggplot(aes(x = n_lmks, y = value,
             linetype = as.factor(ext_mult), color = proc)) +
  facet_grid(stat ~ ., scales = "free_y") +
  geom_line(size = .5) +
  scale_color_manual(values = proc_pal) +
  labs(x = "Number of landmarks", y = "Statistic",
       linetype = "Extension", color = "Procedure") +
  theme(legend.box = "horizontal") ->
  cover_eval_plot_mx
ggsave(here::here("docs/figures/cover-evaluate-mx.pdf"),
       cover_eval_plot_mx,
       width = textwidth * 1/2, height = textwidth / phi)
