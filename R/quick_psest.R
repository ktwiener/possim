

test_settings <- function(prep_pop){
  prep_pop |>
    dplyr::mutate(
      ps_results = purrr::map(data, quick_psest),
      emp_ests = purrr::map(ps_results, emp_est_ps)
    ) |>
    tidyr::unnest(cols = emp_ests) |>
    dplyr::group_by(scenario, effect) |>
    dplyr::summarize(
      across(ends_with("ipt"), ~mean(.x)),
      across(ends_with("smr"), ~mean(.x))
    ) |>
    dplyr::select(scenario, effect, starts_with("ln"), starts_with("r")) |>
    dplyr::mutate(across(starts_with("ln"), ~exp(.x)))
}
quick_psest <- function(dset){

  dset$ps <-
    predict(
      glm(a ~ w,
          data = dset,
          family = binomial),
      newdata = dset,
      type = "response")


  dset %>%
    mutate(
      smr = a*1 + (1-a)*ps/(1-ps),
      ipt = a/ps + (1-a)/(1-ps)
    )

}

emp_est_ps <- function(weights) {
  weights %>%
    dplyr::summarize(
      tot = n(),
      r1 = sum(a*y)/sum(a),
      r0 = sum((1-a)*y)/sum(1-a),
      rr = r1/r0,
      r1ipt = sum(a*y*ipt)/sum(a*ipt),
      r0ipt = sum((1-a)*y*ipt)/sum((1-a)*ipt),
      r1_mipt = (1/tot)*sum(a*y*ipt)/(1/tot*sum(a*ipt)),
      r0_mipt = (1/tot)*sum((1-a)*y*ipt),
      rript = r1ipt/r0ipt,
      r1smr = sum(a*y*smr)/sum(a*smr),
      r0smr = sum((1-a)*y*smr)/sum((1-a)*smr),
      rrsmr = r1smr/r0smr
    )
}

mod_est_ps <- function(pdata){
  browser()
  mods <- list(
  ln_rr_smr = glm(y ~ a, weights = smr, data = pdata, family = binomial(link = "log")),
  ln_or_smr = glm(y ~ a, weights = smr, data = pdata, family = binomial(link = "logit")),

  ln_rr_ipt = glm(y ~ a, weights = ipt, data = pdata, family = binomial(link = "log")),
  ln_or_ipt = glm(y ~ a, weights = ipt, data = pdata, family = binomial(link = "logit"))
  )

  new <- dplyr::tibble(
    a = 0:1
  )

  purrr::imap_dfc(
    mods,
    ~{
      dplyr::tibble(!!.y := exp(.x$coefficients[2]))
    }
  ) %>%
    dplyr::mutate(

    )

}


