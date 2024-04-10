## Results from delicatessen

library(dplyr)
source("R/performance.R")
source("R/utils.R")
source("R/settings.R")

date <- "2024-03-27"
positivity <- c("Full", "Partial")
delt <- c("None","Homogeneous")
probw <- c(25, 50, 75)

measures <- purrr::pmap_dfr(
    .l =  expand.grid(positivity = positivity, delt = delt, probw = probw, date = date, stringsAsFactors = F),
    .f =  function(positivity, delt, probw, date) {

        setting_file <- dplyr::last(list.files(path = "data/settings/", pattern = date, full.names = T))
        patt <- paste0(tolower(c(positivity, delt, probw)), collapse = "-")
        results_files <- list.files(sprintf("data/simulations/%s/", date), pattern = patt, full.names = T)

        results <- read.csv(results_files)
        settings <- readRDS(setting_file)
        effects <- calculate_effects(settings)

        effects <- effects[with(effects, label == paste0(positivity, " exchangeability") & effect == delt & 100*w==probw),]

        hold <- cbind(results[grepl("delta", results$Param), ], effects)

        hold |>
            dplyr::mutate(
            ## Effects
            deltaor = if_else(grepl("smr", Param), att, ate),
            deltarr = deltaor/(1-py0+(py0*deltaor)),
            delta = log(deltarr),

            ## Measures
            cov = LCL <= delta & delta <= UCL,
            mse = (Coef - delta)^2,
            bias = Coef - delta,
            positivity = positivity, 
            efftype = delt, 
            probw = probw

            ) 

    }
)

saveRDS(measures, sprintf("data/results/raw/%s-measures.rds", date))

measures_summary <- measures |>
            dplyr::group_by(positivity, efftype, probw, Param) |>
            dplyr::summarize(
                delta = delta[1],
                expest = mean(Coef),
                bias = mean(bias),
                ese = sqrt(var(Coef)),
                ase = sqrt(mean(Variance)),
                cov = mean(cov),
                mse = mean(mse)
            )

measures_summary  %>% 
    dplyr::transmute(
        trteffect = factor(efftype, levels = c("None", "Homogeneous"), labels = c("No treatment effect", "Homogeneous treatment effect")), 
        scenario = positivity,  
        estimator = if_else(grepl("hajek|smr", Param), "Hajek", "Horvitz-Thompson"),
        weight = if_else(grepl("smr",Param), "SMR", "IPT"),
        pw = probw,     
        across(c(bias, ase, ese, mse, cov), ~formatC(.x, digits = 2, format = "f"))                  
    )  %>% 
    arrange(trteffect, scenario, estimator, weight, pw) %>%
    gt(groupname_col = c("trteffect")) %>%
    cols_label(
        scenario = "Scenario",
        weight = "Weight",
        estimator = "Estimator",
        pw = "P(W=1)",
        bias = "Bias",
        ase = "ASE",
        ese = "ESE",
        mse = "MSE",
        cov = "Coverage"
    ) %>%
  gt::gtsave(paste0(esttype, "-rr-sims.rtf"), "data/results/")

