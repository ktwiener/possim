## Script to run through simulations
library(dplyr)
devtools::load_all()

## Create population settings
pop <- purrr::pmap_dfr(
  .l = settings[1,],
  .f = popgen
)

## Simulation settings
settings_tbl(settings[1,])

## Estimate propensity scores
pspop <- pop |>
  group_by(scenario, sim) |>
  psest()

## Create estimates
estimates <- dplyr::group_modify(.data = pspop,
                 .f = ~cond_effects(.x))

## Calculate variances
unweighted_variance <- estimate_variance(pop)
bootstrap_variance <- estimate_bootstrap_variance(pop, bootstraps = 100)

all_variance <- dplyr::bind_rows(bootstrap_variance, unweighted_variance) |>
  dplyr::arrange(scenario, sim, type)


measures <- performance_measures(estimates, all_variance,
                                 delta = log(1), sims = settings$sims[1])

