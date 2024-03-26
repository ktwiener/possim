# Script to run through simulations ----
setwd("~/Documents/Development/possim")
## Set up environment/libraries ----
library(dplyr)
library(geex)
library (furrr)
library(minpack.lm)
library(stringr)

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
  set_settings,
  effects = c("Homogeneous")
)

settings %>%
  mutate(pa_w1 = pw_a1*pa/as.numeric(pw),
         pa_w0 = (1-pw_a1)*pa/(1-as.numeric(pw)),
         pa_w0_test = (0.2 - 0.25*as.numeric(pw))/(1-as.numeric(pw)),
         pa_w1_test = (0.2 - 0.1*(1-as.numeric(pw)))/as.numeric(pw)
  )


saveRDS(settings,
        file = sprintf("data/settings/%s-sims%s-scen%s-pa%s-settings.rds", Sys.Date(),
                       settings$sims[1], nrow(settings), settings$pa[1]*10))
## ## Simulation settings ----
effects <- calculate_effects(settings)

# settings_tbl(effects) |>
#   write.csv(file = sprintf("results/settings-%s-sims%s-scen%s.csv", Sys.Date(), settings$sims[1], nrow(settings)))

## Create population based on settings ----
set.seed(20130127)

pop <- purrr::pwalk(
  .l = settings[,!names(settings) %in% c("pa1_w0", "pa1_w1")],
  .f = popgen
)



