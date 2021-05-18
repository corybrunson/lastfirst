# session library
library(tidyverse)
library(landmark)

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  mx_data <- "~/Desktop/covid19-mx/data"
  save_dir <- "data/cover"
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  mx_data <- "/blue/rlaubenbacher/jason.brunson/covid19-mx/data"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
} else {
  stop("Cannot recognize working directory.")
}

# numbers of landmarks
ns_lmks <- c(6L, 12L, 24L, 36L, 48L, 60L, 120L)
# multiplicative extensions
exts_mult <- c(0, .1, .2)

# Mexican covid-19 data
mx_cases <- read_rds(file.path(mx_data, "mx-cases.rds"))
# binary predictor and response matrices
mx_pred <- mx_cases %>%
  select(-starts_with("fecha_"), -intubado, -def) %>%
  mutate_all(as.integer) %>%
  as.matrix()
mx_resp <- mx_cases %>%
  select(intubado, def) %>%
  mutate_all(as.integer) %>%
  as.matrix()

for (n_lmks in ns_lmks) for (ext_mult in exts_mult) {
  message("Using ", n_lmks, " landmarks with ", ext_mult, " extension...")
  
  # landmark selection and cover construction
  mx_lmks_mm <- landmarks_maxmin(
    mx_pred, dist_method = "cosine",
    num = n_lmks, seed_index = "minmax",
    cover = TRUE,
    extend_radius = extension(mult = ext_mult)
  )
  mx_lmks_lf <- landmarks_lastfirst(
    mx_pred, dist_method = "cosine",
    num = n_lmks, seed_index = "firstlast",
    cover = TRUE,
    extend_cardinality = extension(mult = ext_mult)
  )
  
  # cover set outcome rates
  mx_lmks_mm %>%
    mutate(
      ordinal = row_number(),
      size = map_int(cover_set, length),
      intubado = map_dbl(cover_set, ~ mean(mx_resp[.x, 1L])),
      def = map_dbl(cover_set, ~ mean(mx_resp[.x, 2L]))
    ) ->
    mx_lmks_mm
  mx_lmks_lf %>%
    mutate(
      ordinal = row_number(),
      size = map_int(cover_set, length),
      intubado = map_dbl(cover_set, ~ mean(mx_resp[.x, 1L])),
      def = map_dbl(cover_set, ~ mean(mx_resp[.x, 2L]))
    ) ->
    mx_lmks_lf
  
  # save covers
  write_rds(mx_lmks_mm, file.path(save_dir, str_c(
    "mx-cover-lmk", n_lmks, "-ext", ext_mult*100, "-mm.rds"
  )))
  write_rds(mx_lmks_lf, file.path(save_dir, str_c(
    "mx-cover-lmk", n_lmks, "-ext", ext_mult*100, "-lf.rds"
  )))
  
}
