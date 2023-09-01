
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
     # smr <- 1 # trt*1 + (1-trt)*pscore/(1-pscore)
     # ef_smr_r1 <- trt*smr*outcome - mu[3]
     # ef_smr_r0 <- (1 - trt)*smr*outcome - mu[4]
     #
     # ef_smr_lnrr <- log(mu[3]) - log(mu[4]) - delta[3]
     # ef_smr_lnor <- log(mu[3]*(1-mu[4])) - log(mu[3]*(1-mu[4])) - delta[4]


    return(c(ests, ef_ipt_r1, ef_ipt_r0, ef_ipt_lnrr, ef_ipt_lnor))#,
             #ef_smr_r1, ef_smr_r0, ef_smr_lnrr, ef_smr_lnor))
  }
}

lm_redo <- function(f, start) {
  x <- minpack.lm::nls.lm(fn = f, par = start,
                          lower = c(-Inf, -Inf, 0, 0, -Inf, -Inf),
                          upper = c(Inf, Inf, 1, 1, Inf, Inf))
  #x$root <- x$par

  x
}

with_log <- geex::m_estimate(
  estFUN = eefun_ipt,
  data = dset,
  root_control = setup_root_control(FUN = lm_redo,
                                    start = c(-2, .5, .5, .5, .5, .5),
                                    roots_name = "par")
  )


