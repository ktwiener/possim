# Simulation settings
#
# This file contains the settings that will be applied to the simulation.

logit <- function(p)log(p/(1-p))
set.seed(20130127)
set_settings <- function(w, effects = c("None", "Homogeneous", "Heterogeneous")){
  tibble::tribble(
    ~label,                    ~effect,         ~delta    ,  ~pw_a1,
    "Full exchangeability"   ,  "None"         ,  log(1)  ,       w,
    "Partial exchangeability",  "None"         ,  log(1)  ,       1,
    "Full exchangeability"   ,  "Homogeneous"  ,  deff,           w,
    "Partial exchangeability",  "Homogeneous"  ,  deff,           1,
    "Full exchangeability"   ,  "Heterogeneous",  11.5,           w,
    "Partial exchangeability",  "Heterogeneous",  11.5,           1
  ) %>%
    dplyr::mutate(
      pw_a0 = (w-pa*pw_a1)/(1-pa),
      sims = sims,
      n = n,
      pa = pa,
      alpha = logit(0.01),
      wbeta = logit(0.08)-alpha,
      #alpha = logit(0.10) - pa*delta - w*wbeta,
      ## I want the ATT to be 0.8 -- log(5) + x = log(0.8) -> x = log(0.8) - log(5)
      eta = if_else(effect != "None", deff - delta, 0),
      pw = as.character(w)
    ) %>%
    dplyr::filter(effect %in% effects)
}


calculate_effects <- function(settings){
  settings |>
    dplyr::mutate(
      w    = as.numeric(pw),
      y0w0 = expit(alpha),
      y1w0 = expit(alpha + delta),
      y0w1 = expit(alpha + wbeta),
      y1w1 = expit(alpha + wbeta + delta + eta),
      py0  = y0w1*pw_a0 + y0w0*(1-pw_a0),
      py1  = y1w1*pw_a1 + y1w0*(1-pw_a1),
      py   = py0*(1-pa) + py1*(pa),
      ate = exp((1-w)*delta + w*(delta+eta)),
      att = exp(pw_a1*(delta+eta) + (1-pw_a1)*delta)
    )
}

settings_tbl <- function(x){
  x |>
    dplyr::select(label, effect, sims, n, pa, alpha, wbeta, delta, eta, pw_a1, pw_a0, y0w0, y1w0, y0w1, y1w1, w0OR, w0RR, w1OR, w1RR) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(cols = c(6:last_col())) |>
    dplyr::mutate(value = dplyr::if_else(grepl("y|delta|OR|RR|alpha", name),
                                         sprintf("%3.2f", value),
                                         as.character(value)),
                  sortn = if_else(grepl("y", name), 5.5, as.numeric(row_number())),
                  name_f = factor(name,
                                levels = setting_lvls,
                                labels = setting_lbls[setting_lvls])
    ) |>
    tidyr::pivot_wider(id_cols = c(effect,name), names_from = label)
}

setting_lbls <- list(
   label = "Scenario label for the simulation settings",
   sims  = "Number of simulations",
   n     = "Sample size",
   pa    = "Prevalence of treatment in the population",
   alpha = "Balancing intercept",
   delta = "Log-odds of treatment effect in W = 0",
   eta   = "Log-odds of change of treatment effect in W = 1 from W = 0",
   pw_a1 = "Prevalence of W = 1 among the treated in the population",
   pw_a0 = "Prevalence of W = 1 among the untreated in the population",
   y0w0  = "$P(Y^{a=0}|W=0)$",
   y0w1  = "$P(Y^{a=0}|W=1)$",
   y1w0  = "$P(Y^{a=1}|W=0)$",
   y1w1  = "$P(Y^{a=1}|W=1)$",
   w1OR = "Treatment effect (OR) among W = 1",
   w1RR = "Treatment effect (RR) among W = 1",
   w0OR = "Treatment effect (OR) among W = 0",
   w0RR = "Treatment effect (RR) among W = 0",
   wbeta = "Log-odds of outcome in those with W = 1"

)

setting_lvls <- c("alpha", "y0w0", "y0w1", "y1w0", "y1w1",
                  "pw_a1", "pw_a0", "delta", "eta", "wbeta", "w1OR", "w0OR", "w1RR", "w0RR")


