# Script to run through simulations ----
setwd("~/Documents/Development/possim")
## Set up environment/libraries ----
library(dplyr)
library(geex)
library (furrr)
library(minpack.lm)
library(stringr)
library(rootSolve)

# Run quick estimation or m-estimation?
quick <- FALSE

# Prevalence of treatment
pa <- 0.2

# Treatment effect
deff <- log(0.8)

# Number of simulations
sims <- 5000
# Number of patients per simulation
n <- 6000
# Effect of W on outcome
#wbeta <- log(3)
# Prevalence of W in treated
pws <- c(0.25, 0.5, 0.75)

all_scripts <- list.files("R", pattern = "*.R", full.names = T)
purrr::walk(all_scripts, source)
#devtools::load_all()

settings <- purrr::map_dfr(
  pws,
  set_settings2,
  effects = c("None", "Homogeneous")
)


saveRDS(settings,
        file = sprintf("data/settings/%s-sims%s-scen%s-pa%s-settings.rds", Sys.Date(),
                       settings$sims[1], nrow(settings), settings$pa[1]*10))

# settings_tbl(effects) |>
#   write.csv(file = sprintf("results/settings-%s-sims%s-scen%s.csv", Sys.Date(), settings$sims[1], nrow(settings)))

## Create population based on settings ----
set.seed(20130127)

pop <- purrr::pwalk(
  .l = settings,
  .f = popgen2
)



