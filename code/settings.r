# constants
phi <- (sqrt(5)+1)/2

# document dimensions ('cm')
textwidth <- 15.11293 * 2

# global options
options(tibble.print_min = 12L, tibble.print_max = 24L)
#options(ggplot2.discrete.colour = "Set1", ggplot2.discrete.fill = "Set1")
ggplot2::theme_set(ggplot2::theme_bw())
ggplot2::theme_update(legend.position = "bottom", legend.box = "vertical")

if (FALSE) {
  `*.unit` <- function (x, y) {
    u <- attr(x, "unit")
    x <- as.numeric(x) * y
    grid::unit(x, u)
  }
  `/.unit` <- function (x, y) {
    u <- attr(x, "unit")
    x <- as.numeric(x) / y
    grid::unit(x, u)
  }
}
