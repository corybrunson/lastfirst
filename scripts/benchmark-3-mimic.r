library(bench)
library(dplyr)
# installed version of landmark package
library(landmark)

# MIMIC-III RT-data sets
if (stringr::str_detect(here::here(), "corybrunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  lastfirst_dir <- here::here()
} else if (stringr::str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  lastfirst_dir <- "~/lastfirst"
} else {
  stop("Cannot recognize working directory.")
}

# R implementations of maxmin and lastfirst using any distance
procedure_benchmark <- function(
  x, dist_method = "euclidean", num = 1L
) {
  mark <- mark(
    maxmin = landmarks_maxmin(
      x, dist_method = dist_method, num = num, engine = "R"
    ),
    lastfirst = landmarks_lastfirst(
      x, dist_method = dist_method, num = num, engine = "R"
    ),
    check = FALSE
  )
  mark <- select(mark, expression:total_time)
  mark <- mutate(
    mark,
    distance = dist_method,
    n = nrow(x),
    procedure = factor(expression, levels = c("maxmin", "lastfirst")),
    implementation = "R",
    num = num
  )
  mark
}

if (file.exists(file.path(lastfirst_dir, "data/mark-mimic.rds"))) {
  readr::read_rds(file.path(lastfirst_dir, "data/mark-mimic.rds")) ->
    lmk_marks
} else {
  # initialize benchmark tests data frame
  bench::mark() %>%
    mutate(expression = as.character(names(expression))) %>%
    select(expression:total_time) ->
    lmk_marks
}

# MIMIC-III care units
careunits <- readr::read_rds(file.path(rt_data, "mimic-units.rds"))
careunits <- setdiff(careunits, c("NICU", "NWARD", unique(lmk_marks$data)))
for (careunit in careunits) {
  message("Benchmarking ", careunit, " data...")
  # binary data table and distance object
  file.path(rt_data, paste0("mimic-", tolower(careunit), "-cases.rds")) %>%
    readr::read_rds() %>%
    select(-subject_id, -hadm_id, -starts_with("mortality")) %>%
    as.matrix() ->
    x
  # procedure benchmarks
  for (denom in c(1024, 512, 256, 128)) {
    # numbers of landmark points
    num <- as.integer(nrow(x) / denom)
    # benchmarks
    marks <- procedure_benchmark(x, dist_method = "cosine", num = num)
    marks <- mutate(marks, data = careunit, distance = "cosine")
    marks <- mutate(marks, expression = as.character(names(expression)))
    print(marks)
    lmk_marks <- bind_rows(lmk_marks, marks)
    readr::write_rds(lmk_marks, file.path(lastfirst_dir, "data/mark-mimic.rds"))
  }
}
