library(tidyverse)
library(tdaunif)
library(landmark)
library(reticulate)
use_condaenv("r-reticulate")
gd <- import("gudhi")

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  save_dir <- "data/cover"
  lastfirst_dir <- here::here()
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
  lastfirst_dir <- "~/lastfirst"
} else {
  stop("Cannot recognize working directory.")
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

# number of iterations of each experiment
n_iter <- 120L
# number of points sampled from trefoil
n_pts <- 240L
# maximum number of landmarks
n_lmk <- 60L

# regular interval cover
interval_cover <- function(x, num) {
  # ranges in each coordinate dimension
  rans <- apply(x, 2L, function(y) range(y))
  # cover set bounds
  
}
# covers to include in comparison
cover_funs <- list(
  
)
