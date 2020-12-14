# session library
library(tidyverse)
library(landmark)

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  save_dir <- "data/cover"
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
} else {
  stop("Cannot recognize working directory.")
}

# numbers of landmarks
ns_lmks <- c(6L, 12L, 24L, 36L, 48L, 60L, 120L)
# multiplicative extensions
exts_mult <- c(0, .1, .2)

# binary data with cross-validation indices
file.path(rt_data, str_c("mimic-tsicu-cases.rds")) %>%
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

for (n_lmks in ns_lmks) for (ext_mult in exts_mult) {
  print(c(n_lmks, ext_mult))
  
  # landmark selection and cover construction
  unit_lmks_mm <- landmarks_maxmin(
    unit_pred, dist_method = "cosine",
    num = n_lmks, seed_index = "minmax",
    cover = TRUE,
    extend_radius = extension(mult = ext_mult)
  )
  unit_lmks_lf <- landmarks_lastfirst(
    unit_pred, dist_method = "cosine",
    num = n_lmks, seed_index = "firstlast",
    cover = TRUE,
    extend_cardinality = extension(mult = ext_mult)
  )
  
  # cover set mortality rates
  unit_lmks_mm %>%
    mutate(
      ordinal = row_number(),
      size = map_int(cover_set, length),
      mort = map_dbl(cover_set, ~ mean(unit_resp[.x, 1L]))
    ) ->
    unit_lmks_mm
  unit_lmks_lf %>%
    mutate(
      ordinal = row_number(),
      size = map_int(cover_set, length),
      mort = map_dbl(cover_set, ~ mean(unit_resp[.x, 1L]))
    ) ->
    unit_lmks_lf
  
  # save covers
  write_rds(unit_lmks_mm, file.path(save_dir, str_c(
    "unit-cover-TSICU-lmk", n_lmks, "-ext", ext_mult*100, "-mm.rds"
  )))
  write_rds(unit_lmks_lf, file.path(save_dir, str_c(
    "unit-cover-TSICU-lmk", n_lmks, "-ext", ext_mult*100, "-lf.rds"
  )))
  
}
