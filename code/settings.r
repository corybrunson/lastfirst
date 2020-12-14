# constants
phi <- (sqrt(5)+1)/2

# document dimensions
textwidth <- grid::unit(15.11293, "cm")

# global options
options(tibble.print_min = 12L, tibble.print_max = 24L)
#options(ggplot2.discrete.colour = "Set1", ggplot2.discrete.fill = "Set1")
ggplot2::theme_set(ggplot2::theme_bw())
ggplot2::theme_update(legend.position = "bottom")
