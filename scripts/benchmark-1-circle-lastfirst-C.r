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

# initialize benchmark tests data frame
bench::mark() %>%
  mutate(expression = as.character(names(expression))) %>%
  select(expression:total_time) ->
  lmk_marks

# lastfirst benchmarks
for (n in c(60L, 360L, 1680L)) {
  for (num in c(6L, 24L, 48L)) {
    set.seed(n)
    x <- tdaunif::sample_circle(n, sd = .2)
    mark <- tryCatch(
      expr = R.utils::withTimeout(bench::mark(
        landmarks_lastfirst(x, num = num, cardinality = NULL, engine = "C++")
      ), timeout = 3600),
      TimeoutException = function(ex) {
        cat("Timeout at `num = ", num, ", n = ", n, "`.\n", sep = "")
      }
    )
    if (is.null(mark)) next
    mark <- select(mark, min:total_time)
    mark <- mutate(
      mark,
      data = "circle",
      distance = "euclidean",
      n = nrow(x),
      procedure = "lastfirst",
      implementation = "C++",
      num = num
    )
    print(mark)
    lmk_marks <- bind_rows(lmk_marks, mark)
    readr::write_rds(
      lmk_marks,
      file.path(lastfirst_dir, "data/mark-circle-lastfirst-C.rds")
    )
  }
}
