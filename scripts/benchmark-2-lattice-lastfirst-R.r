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

# lattice samples
sample_lattice <- function(n, num_x = 6L, num_y = 6L, pmf = function(x, y) 1) {
  m <- outer(seq(num_x) - 1L, seq(num_y) - 1L, pmf)
  prob <- m / sum(m)
  x <- sample(x = num_x * num_y, size = n, prob = prob, replace = TRUE)
  p <- expand.grid(x = seq(num_x) - 1L, y = seq(num_y) - 1L)
  as.matrix(p[x, , drop = FALSE])
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
    x <- sample_lattice(n, num_x = 24L, num_y = 12L,
                        pmf = function(x, y) 2 ^ (- y * y))
    mark <- tryCatch(
      expr = R.utils::withTimeout(bench::mark(
        landmarks_lastfirst(x, num = num, cardinality = NULL, engine = "R")
      ), timeout = 3600),
      TimeoutException = function(ex) {
        cat("Timeout at `num = ", num, ", n = ", n, "`.\n", sep = "")
      }
    )
    if (is.null(mark)) next
    mark <- select(mark, min:total_time)
    mark <- mutate(
      mark,
      data = "lattice",
      distance = "euclidean",
      n = nrow(x),
      procedure = "lastfirst",
      implementation = "R",
      num = num
    )
    print(mark)
    lmk_marks <- bind_rows(lmk_marks, mark)
    readr::write_rds(
      lmk_marks,
      file.path(lastfirst_dir, "data/mark-lattice-lastfirst-R.rds")
    )
  }
}
