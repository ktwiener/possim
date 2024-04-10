## Results
library(gt)
all_scripts <- list.files("R", pattern = "*.R", full.names = T)
purrr::walk(all_scripts, source)
date <- "2024-04-08"

measures <- readRDS(sprintf("data/results/raw/%s-measures.rds", date))

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
