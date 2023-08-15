
#' Create simulated population
#'
#' @description This function creates the full simulated population under a
#' simulation scenario.
#'
#' @param sims number of simulations
#' @param n    sample size for each simulation
#' @param pa   prevalence of treatment in population
#' @param alpha stabilizing intercept
#' @param delta log-odds treatment effect
#' @param pw_a1 prevalence of confounder in the treated
#' @param pw_a0 prevalence of confounder in the untreated
#' @return A simulated population
#' @export
#'
popgen <- function(label, sims, n, pa, alpha, delta, pw_a1, pw_a0){

  ns <- sims * n # Full population size

  dplyr::tibble(
    scenario = label,
    # Patient ID
    rowid = 1:(ns),
    # Simulation number to group by
    sim   = ceiling(rowid/n),
    # Patient number within simulation
    id    = ceiling(rowid/sim),
    # Indicator variable for treatment
    a     = rbern(ns, pa),
    # Indicator for confounder.
    w     = a*rbern(ns, pw_a1) + (1-a)*rbern(ns, pw_a0),
    # Potential outcome for person under no treatment
    y0    = rbern(ns, expit(alpha + w)),
    # Potential outcome for person under treatment
    y1    = rbern(ns, expit(alpha + w + delta)),
    # Observed outcome
    y     = a*y1 + (1-a)*y0
  )

}

#' Bootstrapping population
#'
#' @description This function applies the specified number of bootstrap to the cohort
#' and runs the propensity score model + effect estimates.
#'
#' @param cohort The sample cohort for which to estimate propensity scores
#' @param form   Model formula to estimate propensity scores.
#' @param bootstraps Number of bootstraps
#'
#' @return Bootstrapped variance
#' @export
#'

estimate_bootstrap_variance <- function(cohort, bootstraps = 10){
  nids <- max(cohort$id)
  purrr::map_dfr(1:bootstraps,
                      ~{
                        cohort[cohort$id %in% sample(1:nids, nids, replace = T), ]|>
                          group_by(scenario, sim) |>
                          psest() |>
                          dplyr::group_modify(.data = _,
                                                     .f = ~cond_effects(.x, c("smr", "ipw")))
                        }
                      ) |>
    dplyr::group_by(scenario, sim, type) |>
    dplyr::summarize(
      lnrr_var = var(lnrr),
      lnor_var = var(lnor),
      boots = n()
    )

}

#' Propensity score estimator
#'
#' @description This function uses logistic regression to estimate the
#' propensity score (the probability of receiving treatment based on covariate values.)
#' It also calculates inverse probability of treatment (IPT) and standardized mortality (SMR) weights
#' based on the propensity score.
#'
#' @param cohort The sample cohort for which to estimate propensity scores
#' @param form   Model formula to estimate propensity scores.
#'
#' @return The sample cohort with columns added for the propensity score, stabilized IPT, and SMR weights.
#' @export
#'
psest <- function(cohort, form = a ~ w){

  cohort$ps <- predict(glm(form, data = cohort, family = "binomial"),
                       newdata = cohort,
                       type = "response")

  cohort |>
    dplyr::mutate(
      # Probability of treatment in the cohort
      pa_o = sum(a)/dplyr::n(),
      # Stabilized IPT weights.
      ipw = a*(pa_o/ps) + (1-a)*(1-pa_o)/(1-ps),
      # SMR weights
      smr = a*1 + (1-a)*(ps/(1-ps))
    )
}

#' Empirical treatment effect estimates
#'
#' @description Calculates weighted odds and risk ratios
#'
#' @param cohort The sample cohort with weight columns.
#' @param wght   Which weighting type to use. Options are "smr" and "ipw" and "none"
#'
#' @return tibble with or and rr calculated according to weight specification(s).
#' @export
#'

cond_effects <- function(cohort, wght = c("none", "smr", "ipw")){

  ## Dummy for no weights
  cohort$none <- 1

  purrr::map_dfr(
    wght,
    cond_effect,
    cohort
  )

}

estimate_variance <- function(cohort){
  twoxtwo <- table(cohort$a, cohort$y)
  margins <- rowSums(twoxtwo)
  props <- twoxtwo/rowSums(twoxtwo)

  cohort |>
    dplyr::group_by(scenario, sim, a, y) |>
    dplyr::summarize(
      n = n()
    ) |>
    dplyr::mutate(
      total = sum(n),
      prop  = n/total,
      prop_var = prop[y==0]/(prop[y==1]*total)
    ) |>
    dplyr::group_by(scenario, sim) |>
    dplyr::summarize(
      lnrr_var = sum(prop_var[1]),
      lnor_var = sum(1/n)) |>
    dplyr::mutate(
      type = "none",
      boots = 0
    )

}

#' Empirical treatment effect estimates, one at a time
#'
#' @description Calculates weighted odds and risk ratios according to a specified
#' weight column.
#'
#' @param cohort The sample cohort with weight columns.
#' @param wght   Which weighting type to use.
#'
#' @return tibble with or and rr calculated according to weight specification.
#'
cond_effect <- function(wght, cohort){

  ## Risk for the treated, accounting for weights.
  r1 <- risk_est(y    = cohort$y[cohort$a==1],
                 wght = cohort[[wght]][cohort$a==1])

  ## Risk for the untreated, accounting for weights.
  r0 <- risk_est(y    = cohort$y[cohort$a==0],
                 wght = cohort[[wght]][cohort$a==0])

  effect_tib(wght, r1, r0)
}

#' Smaller helper function bc it's used in two places.
#'
effect_tib <- function(wght, r1, r0){
  dplyr::tibble(
    type = wght,
    lnrr   = log(r1/r0),
    lnor   = log(or_est(r1, r0))
  )

}

risk_est <- function(y, wght) {evt_est(y, wght)/sum(wght)}

evt_est <- function(y, wght) sum(y*wght)

or_est <- function(p1, p0) (p1/(1-p1))/(p0/(1-p0))

#' Outcome model
#'
#' @description This function uses binomial regression to estimate the
#' the treatment effect on the relative scale in the sample cohort.
#'
#' @param cohort The sample cohort for which to estimate propensity scores
#' @param wght   Which weighting type to use. Options are "smr" and "ipw"
#' @param linksy Which link function to use for the binomial regression. Options are
#' "logit" (the default) and "log" (for log-risk regression.)
#'
#' @return The sample cohort with columns added for the propensity score, stabilized IPT, and SMR weights.
#' @export
#'
#'
out_model <- function(cohort, wght, linksy){
  wght_sym <- rlang::ensym(wght)
  broom::tidy(glm(y ~ a, family = binomial(linksy), weights = eval(wght_sym),
                  data = cohort)) |>
    dplyr::mutate(expest = exp(estimate))
}

performance_measures <- function(ests, vars, delta = log(1), sims){

  ests |>
    calculate_bootstrap_ci(vars, delta)|>
    dplyr::group_by(scenario, type) |>
    dplyr::summarize(
      # Measures
      ## 0. Expected value of estimate
      expect_lnrr = mean(lnrr),
      expect_lnor = mean(lnor),
      ## 1. Bias of the estimate
      bias_lnrr   = expect_lnrr-delta,
      bias_lnor   = expect_lnor-delta,
      ## 2. Empirical variance of estimate
      emp_var_lnrr = sum((lnrr - expect_lnrr)^2)/((sims-1)),
      emp_var_lnor = sum((lnor - expect_lnor)^2)/((sims-1)),
      ## 3. MSE
      mse_lnrr = mean((lnrr - delta)^2),
      mse_lnor = mean((lnor - delta)^2),
      ## 4. Coverage
      cov_lnrr = mean(lnrr_cov),
      cov_lnor = mean(lnor_cov),
      ## 5. ASE
      ase_lnrr = sqrt(mean(lnrr_var)),
      ase_lnor = sqrt(mean(lnor_var)),

      # Monte Carlo Variances
      ## 1. MC variance of bias
      mc_var_bias_lnrr = emp_var_lnrr/sims,
      mc_var_bias_lnor = emp_var_lnor/sims,
      ## 2. MC variance of empirical variance
      mc_var_emp_var_lnrr = emp_var_lnrr/(2*(sims-1)),
      mc_var_emp_var_lnor = emp_var_lnor/(2*(sims-1)),
      ## 3. MC variance of MSE
      mc_var_mse_lnrr = sum(((lnrr-delta)^2-mse_lnrr)^2)/(sims*(sims-1)),
      mc_var_mse_lnor = sum(((lnor-delta)^2-mse_lnor)^2)/(sims*(sims-1)),
      ## 4. MC variance of coverage
      mc_var_cov_lnrr = sqrt(cov_lnrr*(1-cov_lnrr)/sims),
      mc_var_cov_lnor = sqrt(cov_lnor*(1-cov_lnor)/sims),
      ## 5. MC variance of average SE
      mc_var_ase_lnrr = sqrt(var(lnrr_var)/(4*sims*ase_lnrr^2)),
      mc_var_ase_lnor = sqrt(var(lnor_var)/(4*sims*ase_lnor^2))
    )


}

calculate_bootstrap_ci <- function(ests, vars, delta = log(1)) {
  ests |>
    dplyr::left_join(vars, by = c("scenario", "sim", "type")) |>
    dplyr::mutate(
      lnrr_lcl = lnrr - 1.96*sqrt(lnrr_var),
      lnrr_ucl = lnrr + 1.96*sqrt(lnrr_var),
      lnor_lcl = lnor - 1.96*sqrt(lnor_var),
      lnor_ucl = lnor + 1.96*sqrt(lnor_var),
      lnrr_cov = lnrr_lcl <= delta & delta <= lnrr_ucl,
      lnor_cov = lnor_lcl <= delta & delta <= lnor_ucl
    ) |>
    ungroup()
}



#' True treatment effect
#'
#' @description Uses the potential outcomes in the cohort to calculate the true
#' effect estimates on the relative scale (odds and risk)
#'
#' @param cohort The sample cohort which has potential outcomes.
#'
#' @return tibble with true OR and RR.
#' @export
#'

true_effects <- function(cohort){
  ## Risk of outcome, if everyone had been treated.
  r1 <- risk_est(cohort$y1, rep(1, nrow(cohort)))

  ## Risk of outcome, if no one had been treated.
  r0 <- risk_est(cohort$y0, rep(1, nrow(cohort)))

  effect_tib("truth", r1, r0)

}
