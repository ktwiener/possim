# Script to create results/output

## Set up environment/libraries ----
library(dplyr)
source("R/performance.R")
source("R/utils.R")
source("R/settings.R")

esttype <- "mest"
effects <- calculate_effects(settings)

## Combine each simulation
exch <- c("full", "partial")
eff <- c("none", "homogeneous")
pos <- c("25", "50", "75")

results_mest <- purrr::pmap_dfr(
  .l = expand.grid(exch, eff, pos, stringsAsFactors = F),
  .f = function(Var1, Var2, Var3){
    filepat <- tolower(sprintf("%s-%s-%s", Var1, Var2, Var3))
    read_chunked(filepat) %>%
      dplyr::mutate(
        scenario = paste0(stringr::str_to_sentence(Var1), " exchangeability"),
        effect = stringr::str_to_sentence(Var2),
        pw = as.character(as.numeric(Var3)/100)
      )
  }
) %>%
  dplyr::full_join(effects %>% mutate(scenario = label),
                   by = c("scenario", "effect", "pw"))

saveRDS(results_mest, sprintf("data/combined/%s-sims%s-scen%s-pa%s-%s.rds", Sys.Date(),
                              settings$sims[1], nrow(settings), settings$pa[1]*10, esttype))




