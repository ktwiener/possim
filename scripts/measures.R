## Getting measures from delicatessen

library(dplyr)
source("R/performance.R")
source("R/utils.R")
source("R/settings.R")

date <- "2024-04-08"
positivity <- c("Full", "Partial")
delt <- c("None","Homogeneous")
probw <- c(25, 50, 75)

measures <-
  purrr::pmap_dfr(
    .l =  expand.grid(positivity = positivity, delt = delt, probw = probw,
                      date = date, stringsAsFactors = FALSE),
    .f =  function(positivity, delt, probw, date) {
      setting_file <- dplyr::last(list.files(path = "data/settings/",
                                             pattern = "2024-04-08", full.names = T))
      patt <- paste0(tolower(c(positivity, delt, probw)), collapse = "-")
      results_files <- list.files(sprintf("data/simulations/%s", date), pattern = patt, full.names = T)
      results <- bind_rows(
        read.csv(results_files[!grepl("crude", results_files)]),
        read.csv(results_files[grepl("crude", results_files)]))

      settings <- readRDS(setting_file)
      effects <- settings %>%
        dplyr::filter(label == paste0(positivity, " exchangeability") & effect == delt & 100*pw==probw) %>%
        dplyr::mutate(py0w1 = expit(yint + wycoef),
                      py0w0 = expit(yint),
                      py0 = pw*py0w1 + (1-pw)*py0w0,
                      deltaor = exp(delta),
                      deltarr = deltaor/(1-py0+(py0*deltaor)))

      hold <- cbind(results[grepl("delta", results$Param), ], effects)

      measures <- hold |>
        dplyr::mutate(
          ## Measures
          cov = LCL <= delta & delta <= UCL,
          mse = (Coef - delta)^2,
          bias = Coef - delta,
          positivity = positivity,
          efftype = delt,
          probw = probw)

      measures
      }
)

saveRDS(measures, sprintf("data/results/raw/%s-measures.rds", date))
