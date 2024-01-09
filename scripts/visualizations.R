## Visualizations
library(ggplot2)
library(wesanderson)
library(dplyr)
library(gt)

all_scripts <- list.files("R", pattern = "*.R", full.names = T)
purrr::walk(all_scripts, source)

esttype <- "mest"

filenm <- dplyr::last(sort(list.files("data/results/raw", pattern = esttype, full.names = T)))

all_measures <-readRDS(filenm)
measures <- all_measures$measures
measures_summary <- performance_summary(measures)



simulation_box(measures, "Homogeneous", T, prefix = esttype) -> p
p +
  theme(
        legend.position = "right"
  )
ggsave(filename = "./data/results/figures/mest-homogeneous-hajek-cirg-plot.png",
       units = "in",
       width = 7.5,
       height = 4
)
simulation_box(measures, "None", T, prefix = esttype)
simulation_box(measures, "Homogeneous", F, prefix = esttype)
ggsave(filename = "./data/results/figures/mest-homogeneous-ht-cirg-plot.png",
       units = "in",
       width = 7.5,
       height = 4
)
simulation_box(measures, "None", F, prefix = esttype)

# simulation_violin(measures, "Homogeneous", T, prefix = esttype)
# simulation_violin(measures, "None", T, prefix = esttype)
# simulation_violin(measures, "Homogeneous", F, prefix = esttype)
# simulation_violin(measures, "None", F, prefix = esttype)

measures_summary %>%
  arrange(desc(effect), scenario, pars, pw) %>%
  mutate(weight = toupper(stringr::str_extract(pars,
                                        pattern = "ipt|smr"))) %>%
  dplyr::filter(!is.na(weight)) %>%
  ungroup %>%
  transmute(
    trteffect = paste0(factor(effect, levels = c("None", "Homogeneous"), labels = c("No treatment effect", "Homogeneous treatment effect")),
                       "(RR = ",  trimws(formatC(exp(delta), digits = 1)), ")"),
    scenario = if_else(lag(scenario)==scenario & !is.na(lag(scenario)), "", scenario),
    weight,
    #weight = if_else(lag(weight)==weight & !is.na(lag(weight)), "", weight),
    estimator = if_else(grepl("hajek|smr", pars), "Hajek", "Horvitz-Thompson"),
    #estimator = if_else(lag(estimator)==estimator & !is.na(lag(estimator)), "", estimator),
    pw,
    across(c(bias, ase, ese, mse, cov), ~formatC(.x, digits = 2, format = "f"))
  ) %>%
  gt(groupname_col = c("trteffect")) %>%
  cols_label(
    scenario = "",
    weight = "",
    estimator = "",
    pw = "P(W=1)",
    bias = "Bias",
    ase = "ASE",
    ese = "ESE",
    mse = "MSE",
    cov = "Coverage"
  ) %>%
  gt::gtsave(paste0(esttype, "-rr-sims.rtf"), "data/results/")

