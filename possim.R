## Script to run through simulations
library(dplyr)
library(geex)
devtools::load_all()

## Create population settings
pop <- purrr::pmap_dfr(
  .l = settings[1,],
  .f = popgen
)

## Simulation settings
settings_tbl(settings) |>
  write.csv(file = "settings.csv")

## Estimate propensity scores
## nest for multiple simulations
ests <- pop |>
  dplyr::filter(sim == 1, scenario == "Full exchangeability") |>
  group_by(scenario, sim) |>
  tidyr::nest(data = c(2, 4:last_col())) |>
  dplyr::mutate(
    ests = purrr::map(data,
                        ~{
                          m <- geex::m_estimate(
                          estFUN = eefun_ipt,
                          data = .x,
                          root_control = setup_root_control(start = c(-2, 0, .5, .5, .5, 0.5, .5, .5)))

                          # dplyr::tibble(
                          #   pars = c("ps_int", "ps_coeff", "risk1", "risk2", "or", "rr"),
                          #   ests = m@estimates,
                          #   vars = diag(m@vcov)
                          # )

                          m
                        }
                      )
    )

ests |>
  select(-data) |>
  ungroup() |>
  tidyr::unnest(cols = ests) |>
  dplyr::group_by(pars)|>
  dplyr::summarize(
    avg_lnest = mean(log(ests)),
    avg_est = exp(avg_lnest)
  )
#
# pspop |>
#   filter(sim == 1) |>
#   distinct(a, w, smr) |>
#   arrange(a,w)
#
# pspop |>
#   filter(sim == 1) |>
#   group_by(a, w) |>
#   tally()

## Create estimates

ests <- pspop |>
  dplyr::mutate(
    estimates = purrr::map(psdata, cond_effects)
  )


## Calculate variances
unweighted_variance <- estimate_variance(pop)
bootstrap_variance  <- estimate_bootstrap_variance(pop, bootstraps = 100)

all_variance <- dplyr::bind_rows(bootstrap_variance) |>
  dplyr::arrange(scenario, sim, type)


measures <- performance_measures(estimates, all_variance,
                                 delta = log(1), sims = settings$sims[1])

measures %>%
  select(type, ends_with("lnrr")) |>
  write.csv("lnrr_measures_partial.csv")

