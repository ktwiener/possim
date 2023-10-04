# Script to create results/output

## Set up environment/libraries ----
library(dplyr)
library(gt)

all_scripts <- list.files("R", pattern = "*.R", full.names = T)
purrr::walk(all_scripts, source)

## Performance measures ----

results <- readRDS(dplyr::last(list.files("data", full.names = T)))

measures <- results |>
  performance_measures() |>
  performance_summary()

measures %>%
  filter(grepl("rr", pars)) %>%
  arrange(desc(delta), scenario) %>%
  ungroup %>%
  transmute(
    trteffect =paste0(rep(c("No treatment effect (RR = ", "Homogeneous treatment effect (RR = "), each = 4),
                      trimws(formatC(exp(delta), digits = 2)), ")"),
    scenario = if_else(lag(scenario)==scenario & !is.na(lag(scenario)), "", scenario),
    weight = toupper(substr(pars, 4, 6)),
    across(c(bias, ase, ese, mse, cov), ~formatC(.x, digits = 2, format = "g"))
  ) %>%
  gt(groupname_col = c("trteffect")) %>%
  cols_label(
    scenario = "",
    weight = "",
    bias = "Bias",
    ase = "ASE",
    ese = "ESE",
    mse = "MSE",
    cov = "Coverage"
  ) %>%
  gt::gtsave("rr_sims.rtf", "data/")
