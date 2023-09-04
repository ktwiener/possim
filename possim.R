# Script to run through simulations ----

## Set up environment/libraries ----
library(dplyr)
library(geex)
library (furrr)
devtools::load_all()
future::plan(multisession)

## ## Simulation settings ----
effects <- calculate_effects(settings)

settings_tbl(effects) |>
  write.csv(file = sprintf("results/settings-%s-sims%s-scen%s.csv", Sys.Date(), settings$sims[1], nrow(settings)))

## Create population based on settings ----
pop <- purrr::pmap_dfr(
  .l = settings,
  .f = popgen
)

## M-estimation for effect estimates and variance.----
### nested for multiple simulations ----
ests <- pop |>
  group_by(scenario, delta, sim) |>
  tidyr::nest(data = c(3, 5:last_col())) |>
  dplyr::mutate(
    ests = purrr::map(data,
                        function(dset) {
                          m <- geex::m_estimate(
                            estFUN = eefun_ipt,
                            data = dset,
                            root_control = setup_root_control(FUN = lm_redo,
                                                              start = c(-2, .5, rep(.5, 8)),
                                                              roots_name = "par")
                          )

                          list(
                            pars = c("ps_beta0", "ps_beta1",
                                     "ef_ipt_r1", "ef_ipt_r0", "ef_ipt_lnrr", "ef_ipt_lnor",
                                     "ef_smr_r1", "ef_smr_r0", "ef_smr_lnrr", "ef_smr_lnor"),
                            ests = m@estimates,
                            vars = diag(m@vcov)
                          )
                        }
                      )
    )

### Join results to settings and save to future proof ----
results <- ests |>
  dplyr::left_join(effects %>% mutate(scenario = label, delta, deltarr = log(deltaRR)),
                   by = c("scenario", "delta")) |>
  select(-data)

saveRDS(results, sprintf("data/%s-sims%s-scen%s-pa%s.rds", Sys.Date(), settings$sims[1], nrow(settings), settings$pa[1]*10))

## Performance measures ----

results |>
  performance_measures() |>
  performance_summary()
