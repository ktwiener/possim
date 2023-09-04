
### Estimating equation for the weighted RR and OR.
eefun_ipt <- function(data){
  confounders <- cbind(1, data$w)
  trt <- data$a
  outcome <- data$y
  # outcome <- data$y
  # trt <- data$a
  # confounder <- data$w

  function(theta) {
    betas <- theta[1:2]
    mu <- theta[3:4]
    delta <- theta[5:6]
    mu_smr <- theta[7:8]
    delta_smr <- theta[9:10]

    # estimate logistic regression parameters
    pscore <-  plogis(confounders %*% betas)
    ests <- (trt - pscore) %*% confounders

    ## IPT weights
    ipt <- trt/pscore + (1-trt)/(1-pscore)
    ef_ipt_r1 <- trt*ipt*outcome - mu[1]
    ef_ipt_r0 <- (1 - trt)*ipt*outcome - mu[2]

    ef_ipt_lnrr <- log(mu[1]) - log(mu[2]) - delta[1]
    ef_ipt_lnor <- log(mu[1]*(1-mu[2])) - log(mu[2]*(1-mu[1])) - delta[2]

    ## SMR weights
    smr <- trt*1 + (1-trt)*pscore/(1-pscore)
    ef_smr_r1 <- trt*(outcome - mu_smr[1])
    ef_smr_r0 <- (1 - trt)*smr*(outcome - mu_smr[2])

    ef_smr_lnrr <- log(mu_smr[1]) - log(mu_smr[2]) - delta_smr[1]
    ef_smr_lnor <- log(mu_smr[1]*(1-mu_smr[2])) - log(mu_smr[2]*(1-mu_smr[1])) - delta_smr[2]


    return(c(ests, ef_ipt_r1, ef_ipt_r0, ef_ipt_lnrr, ef_ipt_lnor,
             ef_smr_r1, ef_smr_r0, ef_smr_lnrr, ef_smr_lnor))
  }
}

lm_redo <- function(f, start) {
  x <- minpack.lm::nls.lm(fn = f, par = start,
                          lower = c(-Inf, -Inf, rep(c(0, 0, -Inf, -Inf), 2)),
                          upper = c( Inf,  Inf, rep(c(1, 1,  Inf,  Inf), 2)),
                          control = nls.lm.control(maxiter = 100))
  #x$root <- x$par

  x
}

# with_ipt <- geex::m_estimate(
#   estFUN = eefun_ipt,
#   data = dset,
#   root_control = setup_root_control(FUN = lm_redo,
#                                     start = c(-2, .5, rep(.5, 8)),
#                                     roots_name = "par")
#   )


