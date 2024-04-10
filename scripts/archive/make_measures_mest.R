

## Performance measures mestimate----

filenm <- dplyr::last(sort(list.files("data/combined", pattern = "-mest", full.names = T)))

results_mest <- readRDS(filenm)

settings <- readRDS(dplyr::last(sort(list.files("data/settings", full.names = T))))

summarize_weights(results_mest, settings)

measures <- results_mest |>
  performance_measures() %>%
  filter(grepl("rr", pars))

measure_summary <- measures |>
  performance_summary()

saveRDS(
  list(measures = measures, measure_summary = measure_summary),
  gsub("-mest", "-mest-measures", gsub("combined", "results/raw", filenm))
)

