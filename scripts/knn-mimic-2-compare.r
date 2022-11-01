# session library
library(tidyverse)
library(tidymodels)

# source and store directories
if (dir.exists("/blue")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  lastfirst_dir <- "~/lastfirst"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
} else if (str_detect(here::here(), "jason.brunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  lastfirst_dir <- here::here()
  save_dir <- "data/cover"
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

# bind all results
file.path(lastfirst_dir, "data") %>%
  list.files(str_c("^auc-(rt|gower|cos)-[a-z]+\\.rds"), full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(measure = str_extract(file, "(rt|gower|cos)")) %>%
  mutate(measure = c(cos = "cosine", gower = "Gower", rt = "RT")[measure]) %>%
  mutate(measure = factor(measure, c("cosine", "Gower", "RT"))) %>%
  mutate(data = map(file, read_rds)) %>%
  select(-file) %>%
  unnest(c(data)) %>%
  # mutate(sampler = fct_inorder(sampler)) %>%
  #mutate(sampler = fct_explicit_na(sampler, na_level = "none")) %>%
  mutate(sampler = factor(
    sampler,
    levels = c("maxmin", "lastfirst", "random")
  )) %>%
  mutate(wt_opt = factor(
    wt_opt,
    levels = c("rank1", "rank2", "triangle", "inverse", "gaussian")
  )) %>%
  print() -> auc_stats
# join care unit sizes
file.path(rt_data) %>%
  list.files("^mimic-[a-zA-Z]+-cases.rds", full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(data = map(file, read_rds)) %>%
  transmute(
    careunit = str_replace(file, "^.*mimic-([a-zA-Z]+)-cases.*$", "\\1"),
    careunit = toupper(careunit),
    size = map_int(data, nrow)
  ) %>%
  right_join(auc_stats, by = "careunit") %>%
  # order care units by size
  mutate(careunit = fct_reorder(careunit, size)) %>%
  print() -> auc_stats

# regress test performance on factors + unit-measure-procedure interactions
auc_stats %>%
  filter(! is.na(sampler)) %>%
  # filter(! is.na(sampler) & measure != "RT") %>%
  transmute(
    AUROC = test_auc, unit = careunit,
    measure = fct_relevel(measure, "RT"),
    sampler = fct_relevel(sampler, "random"),
    landmarks = as.factor(landmarks)
  ) ->
  auc_data
# model with main effects only
auc_data %>%
  lm(formula = AUROC ~ unit + measure + sampler + landmarks) ->
  auc_mod_main
# model with interaction effects
auc_data %>%
  lm(formula = AUROC ~ unit * measure * sampler + landmarks) ->
  auc_mod_prod
# compare main effects
bind_rows(
  mutate(tidy(auc_mod_main), regression = "main effects"),
  mutate(tidy(auc_mod_prod), regression = "interaction effects")
) %>%
  mutate(regression = fct_inorder(regression)) %>%
  filter(term != "(Intercept)" & ! str_detect(term, "\\:")) %>%
  mutate(term = str_remove(term, "unit|measure|sampler")) %>%
  mutate(term = ifelse(
    str_detect(term, "^landmarks[0-9]+$"),
    str_replace(term, "^landmarks([0-9]+)$", "\\1 landmarks"),
    term
  )) %>%
  mutate(term = fct_rev(fct_inorder(term))) %>%
  ggplot(aes(x = estimate, y = term)) +
  # facet_grid(rows = vars(regression)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = estimate - 2 * std.error,
                      xmax = estimate + 2 * std.error,
                      shape = regression)) +
  labs(x = "Main effect", y = NULL, shape = NULL) ->
  auc_reg_main_plot
ggsave(here::here("docs/figures/knn-compare-main.pdf"),
       auc_reg_main_plot,
       width = textwidth, height = textwidth / (phi^2))

# inspect interaction effects
auc_mod_prod %>%
  tidy() %>%
  # filter(str_detect(term, "^unit.+\\:measure.+\\:sample.+$")) %>%
  filter(str_detect(term, "\\:")) %>%
  mutate(term = str_remove_all(term, "unit|measure|sampler")) %>%
  mutate(term = str_replace_all(term, "\\:", " × ")) %>%
  mutate(term = fct_rev(fct_inorder(term))) %>%
  # separate(term, c("unit", "measure", "sampler"), "\\:", remove = FALSE) %>%
  mutate(unit = str_extract(term, "^[A-Z]{3,5}")) %>%
  mutate(measure = str_extract(term, "(cosine|Gower)")) %>%
  mutate(sampler = str_extract(term, "(random|maxmin|lastfirst)$")) %>%
  mutate(interaction = str_count(term, "\\:")) %>%
  ggplot(aes(x = estimate, y = term, shape = measure, color = sampler)) +
  scale_shape(na.value = 1L) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = estimate - 2*std.error,
                      xmax = estimate + 2*std.error)) +
  theme(legend.box = "horizontal") +
  labs(x = "Interaction effect", y = NULL, color = NULL, shape = NULL) ->
  auc_reg_prod_plot
ggsave(here::here("docs/figures/knn-compare-prod.pdf"),
       auc_reg_prod_plot,
       width = textwidth, height = textwidth / phi)
