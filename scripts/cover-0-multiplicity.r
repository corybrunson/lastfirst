# session library
library(tidyverse)
library(landmark)

# source and store directories
if (str_detect(here::here(), "corybrunson")) {
  # laptop
  rt_data <- "~/Desktop/rt-data"
  mx_data <- "~/Desktop/covid19-mx/data"
  lastfirst_dir <- here::here()
  save_dir <- "data/cover"
} else if (str_detect(here::here(), "jason.brunson")) {
  # HiPerGator
  rt_data <- "/blue/rlaubenbacher/jason.brunson/rt-data"
  mx_data <- "/blue/rlaubenbacher/jason.brunson/covid19-mx/data"
  lastfirst_dir <- "~/lastfirst"
  save_dir <- "/blue/rlaubenbacher/jason.brunson/lastfirst/data/cover"
} else {
  stop("Cannot recognize working directory.")
}

# source settings
source(file.path(lastfirst_dir, "code/settings.r"))

# MIMIC-III data

# all care units
careunits <- read_rds(file.path(rt_data, "mimic-units.rds"))
careunits <- setdiff(careunits, c("NICU", "NWARD"))

# initialize aggregate multiplicities data
mult_data <- tibble()

for (careunit in careunits) {
  # binary predictor multiplicities
  file.path(rt_data, str_c("mimic-", tolower(careunit), "-cases.rds")) %>%
    read_rds() ->
    unit_cases
  # table of multiplicities
  unit_cases %>%
    select(-contains("mortality")) %>%
    # subset predictors
    select(
      subject_id, hadm_id,
      starts_with("age"), starts_with("gender"),
      2L + seq(12L)
    ) %>%
    group_by_at(vars(-subject_id, -hadm_id)) %>% count() %>% ungroup() %>%
    select(n) %>%
    group_by(n) %>% count(name = "freq") %>% ungroup() %>%
    mutate(unit = careunit) ->
    unit_mult
  # append to aggregate
  mult_data <- bind_rows(mult_data, unit_mult)
}

# plot unit-bar histograms
mult_data %>%
  mutate(bin = ifelse(n <= 12, as.character(n), "13+")) %>%
  mutate(bin = fct_inseq(bin)) %>%
  group_by(bin, unit) %>% summarize(freq = sum(freq)) %>% ungroup() %>%
  ggplot(aes(x = bin, y = freq)) +
  facet_wrap(~ unit) +
  geom_histogram(stat = "identity", binwidth = 1) +
  labs(x = "Multiplicity", y = "Number of points")

# Mexican COVID-19 data

# load pre-processed data
file.path(mx_data, "covid19mx.rds") %>%
  read_rds() %>%
  # restrict to positive COVID-19 diagnoses
  filter(resultado == "positivo sars-cov-2") %>%
  select(-resultado) %>%
  # death indicator
  mutate(def = ! is.na(fecha_def)) %>%
  # begin after the pandemic hit
  mutate(semana_ingreso = as.integer(lubridate::week(fecha_ingreso))) %>%
  filter(semana_ingreso > 10L) %>%
  arrange(fecha_ingreso) %>%
  select(-starts_with("hasta_"), -starts_with("semana_")) %>%
  # collapse infrequent sectors to 'other'
  mutate(sector = fct_lump_n(sector, n = 60L, other_level = "otra")) %>%
  # replace age with age group
  mutate(edad_a_5 = cut(edad, breaks = seq(-1L, max(edad) + 10L, 10L))) %>%
  select(-edad) %>%
  # drop missing factor levels
  mutate_if(is.factor, fct_drop) %>%
  # exclude some variables
  select(-entidad_res) %>%
  # dummy variables
  fastDummies::dummy_cols(
    c("entidad_um", "sector",
      "sexo", "edad_a_5",
      "tipo_paciente", "nacionalidad"),
    remove_selected_columns = TRUE
  ) %>%
  # logical variables
  mutate_at(vars(intubado, neumonia, diabetes:otro_caso, uci),
            ~ ifelse(. == "si", 1L, 0L)) %>%
  print() -> mx_cases
write_rds(mx_cases, file.path(mx_data, "mx-cases.rds"))

# multiplicity distribution
mx_cases %>%
  select(-starts_with("fecha_"), -intubado, -def) %>%
  group_by_at(vars(everything())) %>% count(name = "mult") %>% ungroup() %>%
  select(mult) %>%
  group_by(mult) %>% count(name = "freq") %>% ungroup() %>%
  print() %>%
  mutate(bin = ifelse(mult <= 12, as.character(mult), "13+")) %>%
  mutate(bin = fct_inseq(bin)) %>%
  group_by(bin) %>% summarize(freq = sum(freq)) %>% ungroup() %>%
  ggplot(aes(x = bin, y = freq)) +
  geom_histogram(stat = "identity", binwidth = 1) +
  labs(x = "Multiplicity", y = "Number of points")
