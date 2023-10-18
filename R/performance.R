# Script to calculate performance measures ----

performance_measures <- function(results) {
  results |>
    dplyr::mutate(
      ests = purrr::map(ests, bind_rows)
    ) |>
    tidyr::unnest(cols = ests) |>
    dplyr::filter(grepl("lnrr|lnor", pars)) |>
    dplyr::mutate(
      deltaor = if_else(grepl("ipt", pars), ate, att),
      deltarr = deltaor/(1-py0+(py0*deltaor)),
      delta = log(if_else(grepl("lnrr", pars), deltarr,deltaor)),
      lnlcl = ests - 1.96*sqrt(vars),
      lnucl = ests + 1.96*sqrt(vars),
      cov = lnlcl <= delta & delta <= lnucl,
      sq.err = (ests - delta)^2,
      std.err = sqrt(vars)
    )
}

performance_summary <- function(x) {
  x |>
    dplyr::group_by(scenario, effect, pars, sims, pw) |>
    dplyr::summarize(
      cov = mean(cov),
      mse = mean(sq.err),
      ase = mean(std.err),
      ese = sd(ests),
      eval = mean(ests),
      expeval = exp(eval)
    ) |>
    dplyr::mutate(bias = eval - delta,
                  mc_se_bias = ese/sqrt(sims))
}
