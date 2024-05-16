## Getting measures from delicatessen

library(dplyr)
source("R/performance.R")
source("R/utils.R")
source("R/settings.R")

date <- "2024-05-02"
positivity <- c("Full", "Partial")
delt <- c("None","Homogeneous")
probw <- c(25, 50, 75)

measures <-
  purrr::pmap_dfr(
    .l =  expand.grid(positivity = positivity, delt = delt, probw = probw,
                      date = date, stringsAsFactors = FALSE),
    .f =  function(positivity, delt, probw, date) {
      setting_file <- dplyr::last(list.files(path = "data/settings/",
                                             pattern = "2024-05-02", full.names = T))
      patt <- paste0(tolower(c(positivity, delt, probw)), collapse = "-")
      results_files <- list.files(sprintf("data/simulations/%s", date), pattern = patt, full.names = T)
      results <- purrr::map_dfr(results_files, readRDS)
        readRDS(results_files[!grepl("crude", results_files)])

      settings <- readRDS(setting_file)
      effects <- settings %>%
        dplyr::filter(label == paste0(positivity, " exchangeability") & effect == delt & 100*pw==probw) %>%
        dplyr::mutate(py0w1 = expit(yint + wycoef),
                      py0w0 = expit(yint),
                      py0 = pw*py0w1 + (1-pw)*py0w0,
                      py1w1 = expit(yint + wycoef + delta),
                      py1w0 = expit(yint + delta),
                      py0 = pw*py0w1 + (1-pw)*py0w0,
                      testpy0 = expit(yint + pw*wycoef),
                      testpy1 = expit(yint + pw*wycoef + delta),
                      py1 = pw*py1w1 + (1-pw)*py1w0,
                      testor = py1*(1-py0)/(py0*(1-py1)),
                      testor2 = testpy1*(1-testpy0)/(testpy0*(1-testpy1)),
                      testrr = py1/py0,
                      testrr2 = testpy1/testpy0,
                      deltaor = exp(delta),
                      deltarr = deltaor/(1-py0+(py0*deltaor)))

      hold <- cbind(results[grepl("lnrr", results$params), ], effects)

      measures <- hold |>
        dplyr::mutate(
          ## Confirm delta
          delta = log(testrr),
          ## Measures
          Param = params,
          Coef = roots,
          LCL = roots - 1.96*se,
          UCL = roots + 1.96*se,
          cov = LCL <= delta & delta <= UCL,
          mse = (roots - delta)^2,
          bias = roots - delta,
          positivity = positivity,
          efftype = delt,
          probw = probw)


      measures
      }
)

saveRDS(measures, sprintf("data/results/raw/%s-measures.rds", date))
