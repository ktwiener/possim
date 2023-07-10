# Simulation settings
#
# This file contains the settings that will be applied to the simulation.

settings <- tibble::tribble(
  ~label,                 ~sims,  ~n, ~pa, ~alpha, ~delta, ~pw_a1, ~pw_a0,
  "Full exchangeability",    10, 100, 0.2,      1, log(1),    0.5,    0.5,
  "Partial exchangeability", 10, 100, 0.2,      1, log(1),      1,    0.5
)

setting_lbls <- list(
   label = "Scenario label for the simulation settings",
   sims  = "Number of simulations",
   n     = "Sample size",
   pa    = "Prevalence of treatment in the population",
   alpha = "Stabilizing intercept to keep probabilities under 1",
   delta = "Treatment effect in the population",
   pw_a1 = "Prevalence of confounder w among the treated in the population",
   pw_a0 = "Prevalence of confounder w among the untreated in the population",
   y0w0  = "$P(Y^{a=0}|W=0)$",
   y0w1  = "$P(Y^{a=0}|W=1)$",
   y1w0  = "$P(Y^{a=1}|W=0)$",
   y1w1  = "$P(Y^{a=1}|W=1)$"
)

settings_tbl <- function(x){

  x |>
    dplyr::mutate(
      y0w0 = expit(alpha),
      y0w1 = expit(alpha + 1),
      y1w0 = expit(alpha + delta),
      y1w1 = expit(alpha + 1 + delta)
    ) |>
    dplyr::group_by(label) |>
    tidyr::pivot_longer(cols = 2:12) |>
    dplyr::mutate(value = dplyr::if_else(grepl("delta|y", name),
                                         sprintf("%3.2f", value),
                                         as.character(value)),
                  sortn = if_else(grepl("y", name), 5.5, as.numeric(row_number())),
                  name = factor(name, levels = names(setting_lbls), labels = setting_lbls)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(sortn, name) |>
    dplyr::select(-sortn, -label) |>
    kableExtra::kable(format = "latex", escape = FALSE,col.names = c(" ", " "), booktabs = T, align = c("lr")) |>
    kableExtra::kable_styling(position = "center", latex_options = "HOLD_position")
}
