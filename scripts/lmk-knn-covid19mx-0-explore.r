
# session library
library(tidyverse)

# load pre-processed data
file.path("~/Desktop/covid19-mx/data/covid19mx.rds") %>%
  read_rds() %>%
  print() -> covid19mx_cases

# data are not sorted by significant dates
covid19mx_cases %>%
  mutate(id = row_number()) %>%
  ggplot(aes(x = id, y = fecha_ingreso)) +
  #ggplot(aes(x = id, y = fecha_sintomas)) +
  geom_point(alpha = .2) +
  stat_smooth()

# symtoms are back-dated from admission
covid19mx_cases %>%
  ggplot(aes(x = fecha_ingreso, y = fecha_sintomas)) +
  geom_point(aes(color = resultado), size = 1.5, shape = 16, alpha = .2)
# deaths are (usually) dated after admission and surge in mid-March
covid19mx_cases %>%
  ggplot(aes(x = fecha_ingreso, y = fecha_def)) +
  geom_point(aes(color = resultado), size = 1.5, shape = 16, alpha = .2)

# testing begins to escalate in mid-March
covid19mx_cases %>%
  mutate(semana_ingreso = lubridate::week(fecha_ingreso)) %>%
  ggplot(aes(x = semana_ingreso, y = after_stat(count), fill = resultado)) +
  geom_bar()
