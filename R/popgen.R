
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
popgen <- function(label, effect, sims, n, pa, alpha, delta, eta, pw_a1, pw_a0, wbeta){

  ns <- sims * n # Full population size

  dplyr::tibble(
    scenario = label,
    effect = effect,
    delta = delta,
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
    y0    = rbern(ns, expit(alpha + wbeta*w)),
    # Potential outcome for person under treatment
    y1    = rbern(ns, expit(alpha + wbeta*w + delta + w*eta)),
    # Observed outcome
    y     = a*y1 + (1-a)*y0
  )

}
