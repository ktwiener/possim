library(dplyr)
library(geex)
library (furrr)
library(minpack.lm)
library(stringr)

source("R/each_pop.R")
source("R/mestimate.R")
source("R/utils.R")
## M-estimation for effect estimates and variance.----
### nested for multiple simulations ----

date <- dplyr::last(list.files("data/population/"))
popnms <- list.files(paste0("data/population/", date), full.names = T)

big_loc <- paste0("data/simulations/", date)
sim_loc <- paste0(big_loc, "/parts")
dir.create(big_loc)
dir.create(sim_loc)

purrr::walk(
  popnms,
  each_pop,
  sim_number = 11:200,
  done_files = list.files(sim_loc),
  sim_loc = sim_loc
)

# ### Join results to settings and save to future proof ----
  # results <- prep_pop |>
  #   dplyr::left_join(effects %>% mutate(scenario = label),
  #                    by = c("scenario", "effect", "pw")) |>
  #   select(-data)
  #
  # saveRDS(results, sprintf("data/%s-sims%s-scen%s-pa%s.rds", Sys.Date(), settings$sims[1], nrow(settings), settings$pa[1]*10))



# future::plan(multisession)
# start <- Sys.time()
# sim_results <- test_settings(prep_pop)
#
# saveRDS(sim_results, sprintf("data/%s-sims%s-scen%s-pa%s-quick.rds", Sys.Date(), settings$sims[1], nrow(settings), settings$pa[1]*10))
#
# ests <- sum_results(sim_results$sim_ests)
# rr_ests <- sum_rr(ests)
# furrr_test_time <- Sys.time() - start
# check$emp_ests
