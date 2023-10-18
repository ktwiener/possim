# Script to create results/output

## Set up environment/libraries ----
library(dplyr)
library(gt)
library(ggplot2)
library(wesanderson)

all_scripts <- list.files("R", pattern = "*.R", full.names = T)
purrr::walk(all_scripts, source)

## Performance measures ----

results <- readRDS(dplyr::last(list.files("data", full.names = T)))

measures <- results |>
  performance_measures() %>%
  filter(grepl("rr", pars))

measures %>%
  select(scenario, effect, pw, sims, delta, pars, ests) %>%
  filter(effect == "None", !grepl("hajek", pars)) %>%
  dplyr::mutate(pars = factor(pars,
                              levels = c("ef_ipt_hajek_lnrr", "ef_ipt_lnrr", "ef_smr_lnrr"),
                              labels = c("IPT (Horvitz Thompson)", "IPT (Hajek)", "SMR"))) %>%
  dplyr::mutate(pw = paste0("P(W=1) = ",pw)) %>%
  ggplot(mapping = aes(y = ests, x = pars, color = pars)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  geom_point(position = "jitter", alpha = 0.5) +
  facet_grid(vars(scenario), vars(pw), switch = "x") +
  scale_color_manual(values=wes_palette(n=2, name="Royal1")) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA)) +
  ylab("log RR") +
  labs(color = "Estimator")

ggsave("data/noneplot.png")

measures %>%
  select(scenario, effect, pw, sims, delta, pars, ests) %>%
  filter(effect == "Homogeneous", !grepl("hajek", pars)) %>%
  dplyr::mutate(pars = factor(pars,
                              levels = c("ef_ipt_hajek_lnrr", "ef_ipt_lnrr", "ef_smr_lnrr"),
                              labels = c("IPT (Horvitz Thompson)", "IPT (Hajek)", "SMR"))) %>%
  dplyr::mutate(pw = paste0("P(W=1) = ",pw)) %>%
  ggplot(mapping = aes(y = ests, x = pars, color = pars)) +
  geom_boxplot() +
  geom_hline(yintercept = log(0.8)) +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  geom_point(position = "jitter", alpha = 0.5) +
  facet_grid(vars(scenario), vars(pw), switch = "x") +
  scale_color_manual(values=wes_palette(n=2, name="Royal1")) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA)) +
  ylab("log RR") +
  labs(color = "Estimator")

ggsave("data/homogeneousplot.png")

measure_summary <- measures |>
  performance_summary()

measure_summary %>%
  arrange(desc(effect), scenario, pars, pw) %>%
  ungroup %>%
  transmute(
    trteffect = paste0(factor(effect, levels = c("None", "Homogeneous"), labels = c("No treatment effect", "Homogeneous treatment effect")),
                       "(RR = ",  trimws(formatC(exp(delta), digits = 1)), ")"),
    scenario = if_else(lag(scenario)==scenario & !is.na(lag(scenario)), "", scenario),
    weight = toupper(substr(pars, 4, 6)),
    #weight = if_else(lag(weight)==weight & !is.na(lag(weight)), "", weight),
    estimator = if_else(grepl("hajek|smr", pars), "Hajek", "Horvitz-Thompson"),
    #estimator = if_else(lag(estimator)==estimator & !is.na(lag(estimator)), "", estimator),
    pw,
    across(c(bias, ase, ese, mse, cov), ~formatC(.x, digits = 2, format = "g"))
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
  gt::gtsave("rr_sims.rtf", "data/")
