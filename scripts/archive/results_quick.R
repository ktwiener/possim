# Script to create results/output

## Set up environment/libraries ----
library(dplyr)
source("R/performance.R")
source("R/utils.R")
source("R/settings.R")

esttype <- "quick"
effects <- calculate_effects(settings)

results_quick <- readRDS(dplyr::last(list.files(path = "./data/", pattern = "quick", full.names = T)))

hold <- results_quick$sim_ests %>%
  left_join(effects, by = c("scenario" = "label", "effect", "pw"))

## Performance measures mestimate----

measures <- hold |>
  performance_measures_quick() %>%
  filter(grepl("rr", name)) %>%
  dplyr::rename(pars = name, ests = value)

measure_summary <- measures |>
  performance_summary()

source("R/visualizations.R")

