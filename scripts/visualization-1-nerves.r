# session library
library(tidyverse)
#library(tidymodels)
library(landmark)
#remotes::install_github("peekxc/simplextree",
#                        ref = "6e34926b3e6c7c990226e23ba075e29e9a374365")
library(simplextree)
# source directories
rt_data <- "~/Desktop/rt-data"
# all care units
careunits <- read_rds(file.path(rt_data, "mimic-units.rds"))
careunits <- setdiff(careunits, c("NICU", "NWARD"))

careunit <- "CCU"
n_lmks <- 24L

# function to construct a nerve from a cover
nerve <- function(cover) {
  # initialize SC
  st <- simplex_tree()
  # add a vertex for each cover set
  for (s in seq_along(cover)) st$insert(s)
  # add a higher-dimensional simplex for each overlap
  ids <- sort(unique(unlist(cover)))
  pb <- progress::progress_bar$new(total = length(ids))
  for (id in ids) {
    pb$tick()
    simp <- which(sapply(cover, `%in%`, x = id))
    if (length(simp) < 2L) next
    if (! st$find(s)) st$insert(simp)
  }
}

# binary data with cross-validation indices
file.path(rt_data, str_c("mimic-", tolower(careunit), "-cases.rds")) %>%
  read_rds() %>%
  print() -> unit_cases
# binary predictor and response matrices
unit_pred <- unit_cases %>%
  select(-subject_id, -hadm_id, -contains("mortality")) %>%
  mutate_all(as.integer) %>%
  as.matrix()
unit_resp <- unit_cases %>%
  select(mortality_hosp) %>%
  mutate_all(as.integer) %>%
  as.matrix()
# landmark selection and cover construction
unit_lmks_mm <- landmarks_maxmin(unit_pred, dist_method = "cosine",
                                 num = n_lmks, seed_index = "minmax",
                                 cover = TRUE,
                                 extend_radius = extension(mult = .25))
unit_lmks_lf <- landmarks_lastfirst(unit_pred, dist_method = "cosine",
                                    num = n_lmks, seed_index = "firstlast",
                                    cover = TRUE,
                                    extend_cardinality = extension(mult = .25))

# visualize and annotate overlapping covers
(unit_nerve_mm <- nerve(unit_lmks_mm))
(unit_nerve_lf <- nerve(unit_lmks_lf))

# save covers and nerves
save(unit_lmks_mm, unit_lmks_lf, unit_nerve_mm, unit_nerve_lf,
     file = here::here(str("data/cover-nerve-", careunit, ".rda")))
