
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

#' Smaller helper function bc it's used in two places.
#'
effect_tib <- function(wght, r1, r0){
  dplyr::tibble(
    type = wght,
    rr   = r1/r0,
    or   = or_est(r1, r0)
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


