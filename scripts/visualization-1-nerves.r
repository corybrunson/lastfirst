# session library
library(tidyverse)
#library(tidymodels)
library(landmark)
#remotes::install_github("peekxc/simplextree",
#                        ref = "6e34926b3e6c7c990226e23ba075e29e9a374365")
library(simplextree)
# source directory
#rt_data <- "~/Desktop/rt-data" # laptop
rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data" # HPG
# store directory
save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst"

# all care units
careunits <- read_rds(file.path(rt_data, "mimic-units.rds"))
careunits <- setdiff(careunits, c("NICU", "NWARD"))
# numbers of landmarks
ns_lmks <- c(12L, 24L, 60L)
# multiplicative extensions
exts_mult <- c(0, .1, .25)

if (FALSE) {
  # custom function to construct a nerve from a cover
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
}

for (careunit in careunits) {
  
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
  
  for (n_lmks in ns_lmks) for (ext_mult in exts_mult) {
    print(c(careunit, n_lmks, exts_mult))
    
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
    
    # visualize and annotate overlapping covers
    (unit_nerve_mm <- nerve(simplex_tree(), unit_lmks_mm$cover, k = 4L))
    (unit_nerve_lf <- nerve(simplex_tree(), unit_lmks_lf$cover, k = 4L))
    
    # save covers and nerves
    write_rds(unit_lmks_mm, file.path(save_dir, "data", str_c(
      "unit-cover-", careunit, "-lmk", n_lmks, "-ext", ext_mult*100, "-mm.rds"
    )))
    write_rds(unit_lmks_lf, file.path(save_dir, "data", str_c(
      "unit-cover-", careunit, "-lmk", n_lmks, "-ext", ext_mult*100, "-lf.rds"
    )))
    write_rds(as.list(maximal(unit_nerve_mm)), file.path(save_dir, "data", str_c(
      "unit-nerve-", careunit, "-lmk", n_lmks, "-ext", ext_mult*100, "-mm.rds"
    )))
    write_rds(as.list(maximal(unit_nerve_lf)), file.path(save_dir, "data", str_c(
      "unit-nerve-", careunit, "-lmk", n_lmks, "-ext", ext_mult*100, "-lf.rds"
    )))
    
  }
  
}
