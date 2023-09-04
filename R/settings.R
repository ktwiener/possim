# Simulation settings
#
# This file contains the settings that will be applied to the simulation.

settings <- tibble::tribble(
  ~label,                 ~sims,  ~n, ~pa, ~alpha, ~delta,  ~pw_a1, ~pw_a0,
  "Full exchangeability",    10, 6000, 0.2,      -2, log(1),   0.5,    0.5,
  "Partial exchangeability", 10, 6000, 0.2,      -2, log(1),     1,    0.5,
  "Full exchangeability",    10, 6000, 0.2,      -2, log(0.8), 0.5,    0.5,
  "Partial exchangeability", 10, 6000, 0.2,      -2, log(0.8),   1,    0.5
)

setting_lbls <- list(
   label = "Scenario label for the simulation settings",
   sims  = "Number of simulations",
   n     = "Sample size",
   pa    = "Prevalence of treatment in the population",
   alpha = "Balancing intercept",
   delta = "Treatment effect in the population (ratio scale)",
   pw_a1 = "Prevalence of W = 1 among the treated in the population",
   pw_a0 = "Prevalence of W = 1 among the untreated in the population",
   y0w0  = "$P(Y^{a=0}|W=0)$",
   y0w1  = "$P(Y^{a=0}|W=1)$",
   y1w0  = "$P(Y^{a=1}|W=0)$",
   y1w1  = "$P(Y^{a=1}|W=1)$",
   deltaOR = "Treatment effect (OR)",
   deltaRR = "Treatment effect (RR)"
)

setting_lvls <- c("alpha", "y0w0", "y0w1", "y1w0", "y1w1",
                  "pw_a1", "pw_a0", "deltaOR", "deltaRR")

calculate_effects <- function(x){
  x |>
    dplyr::mutate(
      y0w0 = expit(alpha),
      y0w1 = expit(alpha + 0.5),
      y1w0 = expit(alpha + delta),
      y1w1 = expit(alpha + 0.5 + delta),
      deltaOR = exp(delta),
      deltaRR = y1w1/y0w1
    )
}

settings_tbl <- function(x){
 x |>
    dplyr::select(alpha, delta, c(9:last_col())) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(cols = c(1, 3:8)) |>
    dplyr::mutate(value = dplyr::if_else(grepl("y|delta", name),
                                         sprintf("%3.2f", value),
                                         as.character(value)),
                  sortn = if_else(grepl("y", name), 5.5, as.numeric(row_number())),
                  name = factor(name,
                                levels = setting_lvls,
                                labels = setting_lbls[setting_lvls])
    ) |>
    tidyr::pivot_wider(id_cols = name, names_from = delta)
}
