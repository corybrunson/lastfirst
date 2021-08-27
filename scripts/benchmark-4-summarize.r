library(tidyverse)

# source and store directories
if (stringr::str_detect(here::here(), "corybrunson")) {
  # laptop
  machine <- "Cory's MacBook Air"
  rt_data <- "~/Desktop/rt-data"
  lastfirst_dir <- here::here()
  library(landmark)
} else if (stringr::str_detect(here::here(), "Users/jason.brunson")) {
  # desktop
  machine <- "Cory's UF iMac"
  rt_data <- "~/Desktop/rt-data"
  lastfirst_dir <- here::here()
  devtools::load_all("~/Documents/proj-active/tda/landmark/")
} else if (stringr::str_detect(here::here(), "home/jason.brunson")) {
  # HiPerGator
  machine <- "HiPerGator cluster"
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  lastfirst_dir <- "~/lastfirst"
  library(landmark)
} else {
  stop("Cannot recognize working directory.")
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

# circle samples
file.path(lastfirst_dir, "data") %>%
  list.files("^mark-circle", full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(data = map(file, read_rds)) %>%
  mutate(file = basename(file)) %>%
  unnest(c(data)) %>%
  select(-expression) %>%
  print() -> mark_circle
# lattice samples
file.path(lastfirst_dir, "data") %>%
  list.files("^mark-lattice", full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  mutate(data = map(file, read_rds)) %>%
  mutate(file = basename(file)) %>%
  unnest(c(data)) %>%
  select(-expression) %>%
  print() -> mark_lattice
# MIMIC-III care units
mark_mimic <- read_rds(file.path(lastfirst_dir, "data/mark-mimic.rds"))

# evaluate circle and lattice experiments
bind_rows(mark_circle, mark_lattice) %>%
  mutate_at(vars(min, median, mem_alloc), as.numeric) %>%
  mutate(
    data = factor(data, levels = c("circle", "lattice")),
    procedure = factor(procedure, levels = c("maxmin", "lastfirst")),
    implementation = factor(implementation, levels = c("original", "C++", "R"))
  ) %>%
  pivot_longer(c(median, mem_alloc), names_to = "stat", values_to = "value") %>%
  mutate(data = fct_inorder(data)) %>%
  mutate(data = fct_recode(data, Circle = "circle", Lattice = "lattice")) %>%
  mutate(stat = fct_inorder(stat)) %>%
  mutate(stat = fct_recode(stat,
                           `Median execution time` = "median",
                           `Total memory allocation` = "mem_alloc")) %>%
  select(data, implementation, n, num, procedure, stat, value) %>%
  print() -> benchmark_circle_lattice
benchmark_circle_lattice %>%
  filter(implementation != "original") %>%
  pivot_wider(names_from = "procedure", values_from = "value") %>%
  mutate(factor = lastfirst / maxmin) %>%
  group_by(implementation) %>%
  summarize(factor = median(factor, na.rm = TRUE))
benchmark_circle_lattice %>%
  ggplot(aes(x = n, y = value,
             group = interaction(procedure, num, implementation))) +
  facet_grid(stat ~ data, scales = "free_y") +
  geom_line(aes(color = implementation, linetype = procedure)) +
  geom_point(aes(color = implementation, shape = procedure)) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Size of point cloud", y = "Benchmark statistic",
       color = "Implementation",
       linetype = "Procedure", shape = "Procedure") ->
  benchmark_circle_lattice_plot
ggsave(here::here("docs/figures/benchmark-circle-lattice.pdf"),
       benchmark_circle_lattice_plot,
       width = textwidth * 2/3, height = textwidth / phi, units = "cm")

# evaluate MIMIC-III experiments
mark_mimic %>%
  mutate_at(vars(min, median, mem_alloc), as.numeric) %>%
  mutate_at(vars(procedure, data), fct_inorder) %>%
  pivot_longer(c(median, mem_alloc), names_to = "stat", values_to = "value") %>%
  mutate(stat = fct_inorder(stat)) %>%
  mutate(stat = fct_recode(stat,
                           `Median execution time` = "median",
                           `Total memory allocation` = "mem_alloc")) %>%
  select(data, n, num, procedure, stat, value) %>%
  print() -> benchmark_mimic
benchmark_mimic %>%
  group_by(stat) %>%
  mutate(max_value = max(value, na.rm = TRUE),
         min_value = min(value, na.rm = TRUE),
         label_value = max_value * exp(.15 * log(max_value / min_value))) %>%
  group_by(data, stat) %>%
  mutate(label = ifelse(row_number() == 1L, as.character(data), NA)) %>%
  ungroup() %>%
  # rank of `num`
  group_by(data, procedure, stat) %>%
  mutate(num_rank = as.factor(rank(num))) %>%
  ggplot(aes(x = n, y = value,
             group = interaction(procedure, num_rank))) +
  facet_grid(stat ~ ., scales = "free_y") +
  geom_line(aes(color = num_rank, linetype = procedure)) +
  geom_point(aes(color = num_rank, shape = procedure)) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set2") +
  geom_text(aes(x = n, y = label_value, label = label),
            size = 3, angle = 90, hjust = 1) +
  labs(x = "Size of cohort", y = "Benchmark statistic",
       linetype = "Procedure", shape = "Procedure") ->
  benchmark_mimic_plot
ggsave(here::here("docs/figures/benchmark-mimic.pdf"),
       benchmark_mimic_plot,
       width = textwidth * 1/3, height = textwidth / phi, units = "cm")
