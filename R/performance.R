# Script to calculate performance measures ----

performance_measures <- function(results) {
  results |>
    dplyr::mutate(
      ests = purrr::map(ests, bind_rows)
    ) |>
    tidyr::unnest(cols = ests) |>
    dplyr::filter(grepl("lnrr|lnor", pars)) |>
    dplyr::mutate(
      lnlcl = ests - 1.96*sqrt(vars),
      lnucl = ests + 1.96*sqrt(vars),
      delta = if_else(grepl("lnrr", pars), deltarr, delta),
      cov = lnlcl <= delta & delta <= lnucl,
      sq.err = (ests - delta)^2,
      std.err = sqrt(vars)
    )
}

performance_summary <- function(x) {
  x |>
    dplyr::group_by(scenario, delta, pars, sims) |>
    dplyr::summarize(
      cov = mean(cov),
      mse = mean(sq.err),
      ase = mean(std.err),
      ese = sd(ests),
      eval = mean(ests)
    ) |>
    dplyr::mutate(bias = eval - delta,
                  mc_var_bias = (ese^2)/sims)
}
