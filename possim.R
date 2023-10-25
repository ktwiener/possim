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
sims <- 200
# Number of patients per simulation
n <- 6000
# Effect of W on outcome
#wbeta <- log(3)
# Prevalence of W in treated
pws <- c(0.25, 0.5, 0.75)

all_scripts <- list.files("R", pattern = "*.R", full.names = T)
purrr::walk(all_scripts[-5], source)
#devtools::load_all()

settings <- purrr::map_dfr(
  pws,
  set_settings,
  effects = c("None", "Homogeneous")
)

## ## Simulation settings ----
effects <- calculate_effects(settings)

# settings_tbl(effects) |>
#   write.csv(file = sprintf("results/settings-%s-sims%s-scen%s.csv", Sys.Date(), settings$sims[1], nrow(settings)))

## Create population based on settings ----

pop <- purrr::pmap_dfr(
  .l = settings,
  .f = popgen
)


## M-estimation for effect estimates and variance.----
### nested for multiple simulations ----
prep_pop <- pop |>
  group_by(scenario, effect, pw, sim) |>
  tidyr::nest() |>
  dplyr::mutate(
    filename = tolower(sprintf("%s-%s-%s-%s", word(scenario, 1), effect, as.numeric(pw)*100, sim))
  )

if (quick){
  future::plan(multisession)
  start <- Sys.time()
  sim_results <- test_settings(prep_pop)

  saveRDS(sim_results, sprintf("data/%s-sims%s-scen%s-pa%s-quick.rds", Sys.Date(), settings$sims[1], nrow(settings), settings$pa[1]*10))

  ests <- sum_results(sim_results$sim_ests)
  rr_ests <- sum_rr(ests)
  furrr_test_time <- Sys.time() - start
  check$emp_ests
} else {
  future::plan(multisession)
  start <- Sys.time()
  furrr::future_walk2(prep_pop$data,
                        prep_pop$filename,
                        function(dset, filenm) {
                          m <- geex::m_estimate(
                            estFUN = eefun_ipt,
                            data = as.data.frame(dset),
                            root_control = setup_root_control(FUN = lm_redo,
                                                          start = c(-2, .5, rep(.5, 4*3)),
                                                          roots_name = "par"))
                          ret <- list(
                            pars = c("ps_beta0", "ps_beta1",
                                     "ef_ipt_r1", "ef_ipt_r0", "ef_ipt_lnrr", "ef_ipt_lnor",
                                     "ef_ipt_hajek_r1", "ef_ipt_hajek_r0", "ef_ipt_hajek_lnrr", "ef_ipt_hajek_lnor",
                                     "ef_smr_r1", "ef_smr_r0", "ef_smr_lnrr", "ef_smr_lnor"),
                            ests = m@estimates,
                            vars = diag(m@vcov)
                            )

                          saveRDS(ret, sprintf("data/simulations/%s.rds", filenm))
                    }
  )

  end <- Sys.time()
  furrr_time2 <- end - start

  # ### Join results to settings and save to future proof ----
  # results <- prep_pop |>
  #   dplyr::left_join(effects %>% mutate(scenario = label),
  #                    by = c("scenario", "effect", "pw")) |>
  #   select(-data)
  #
  # saveRDS(results, sprintf("data/%s-sims%s-scen%s-pa%s.rds", Sys.Date(), settings$sims[1], nrow(settings), settings$pa[1]*10))

}
