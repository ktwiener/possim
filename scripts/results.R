## Results
library(gt)
library(ggplot2)
all_scripts <- list.files("R", pattern = "*.R", full.names = T)
purrr::walk(all_scripts, source)
date <- "2024-05-02"

measures <- readRDS(sprintf("data/results/raw/%s-measures.rds", date))

mcse <- measures |>
  dplyr::group_by(positivity, efftype, probw, Param) |>
  dplyr::mutate(expest = mean(Coef), nsim = n()) |>
  dplyr::summarize(
    mcvar_bias = (1/(nsim[1]*(nsim[1]-1)))*sum((Coef-expest)^2)
  ) |>
  dplyr::mutate(
    mcse_bias = sqrt(mcvar_bias)
  ) %>%
  pull(mcse_bias)

print(paste0("Monte-Carlo SE: ", min(mcse), max(mcse)))
## Table
gttbl <- make_table(measures)

gttbl |>
  gt::gtsave(paste0(date, "-rr.rtf"), "data/results/tables/")

## Box Plots

boxplot(measures)


ggsave(filename = sprintf("data/results/figures/%s-homogeneous-hajek-box.png", date),
       units = "in",
       width = 6,
       height = 6
)

boxplot(measures, effect = "None")

ggsave(filename = sprintf("data/results/figures/%s-none-hajek-box.png", date),
       units = "in",
       width = 6,
       height = 6
)

boxplot(measures, ipttype = "ht")

ggsave(filename = sprintf("data/results/figures/%s-homogeneous-ht-box.png", date),
       units = "in",
       width = 6,
       height = 6
)
