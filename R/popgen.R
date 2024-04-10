
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
popgen <- function(label, effect, pw, sims, n, pa, alpha, delta, eta, pw_a1, pw_a0, wbeta){

  ns <- sims * n # Full population size

  dplyr::tibble(
    scenario = label,
    effect = effect,
    pw = pw,
    # Simulation number to group by
    sim   = rep(1:sims, each = n),
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
  ) -> pop

  dir.create(sprintf("data/population/%s", Sys.Date()), showWarnings = F)

  saveRDS(pop, sprintf("data/population/%s/%s-%s-%s.rds", Sys.Date(),
                       tolower(stringr::word(label, 1)), tolower(effect),
                       100*as.numeric(pw)))

}


## Using balancing intercepts

popgen2 <- function(label, effect, pw, sims, n, pa, awcoef, delta, wycoef, yint){

  ns <- sims * n # Full population size

  ## Generate W
  w <- rbern(ns, pw)


  if (awcoef != Inf){
    ## For full exchangeability
    ## Find root for exposure model
    gen_pr_a <- gen_pr_a_wrapper(w, awcoef)
    roota <- multiroot(f=tosolve, start=c(-1), fxname=gen_pr_a, goalmarg=pa, atol=1e-12)$root

    ## Generate A
    a <- rbinom(ns, 1, gen_pr_a(roota))
  } else{
    ## For partial exchangeability.
    pa_w1 = pa/pw
    pa_w0 = 0
    roota <- NA
    a <- rbern(ns, w*pa_w1 + (1-w)*pa_w0)
  }

  pop <- dplyr::tibble(
    scenario = label,
    effect = effect,
    pw = pw,
    delta = delta,
    ## Simulation number to group by
    sim   = rep(1:sims, each = n),
    a = a,
    w = w,
    y0 = rbern(ns, plogis(yint + wycoef*w)),
    y1 = rbern(ns, plogis(yint + wycoef*w + delta)),
    y = a*y1 + (1-a)*y0,
    roota = roota
  )

  dir.create(sprintf("data/population/%s", Sys.Date()), showWarnings = F)

  write.csv(pop, sprintf("data/population/%s/%s-%s-%s.csv", Sys.Date(),
                       tolower(stringr::word(label, 1)), tolower(effect),
                       100*as.numeric(pw)),
            row.names = FALSE)

}

tosolve <- function(ints, fxname, goalmarg){
  values <- fxname(ints)                          # simulated individual expected values at given intercept value

  if(is.null(dim(values))){                        # obtain mean of simulated individual expected values
    empmeans <-  mean(values)
  }else{
    empmeans <-  colMeans(values)
  }

  diff = empmeans - goalmarg                      # difference from empirical mean to goal marginal distribution
  return(diff)
}

