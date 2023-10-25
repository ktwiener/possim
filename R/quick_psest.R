

test_settings <- function(prep_pop){

  prep_pop |>
    dplyr::mutate(
      ps_results = furrr::future_map(data, quick_psest)
      )-> hold

  hold |>
    dplyr::mutate(
      emp_ests = furrr::future_map(ps_results, emp_est_ps)
    ) |>
    tidyr::unnest(cols = emp_ests) |>
    select(-c(data, ps_results, tot)) |>
    tidyr::pivot_longer(cols = c(r1, r0, lnrr,
                                 n1ipt, n0ipt,
                                 r1ipt, r0ipt, lnrr_hajek_ipt,
                                 r1_mipt, r0_mipt, lnrr_ht_ipt,
                                 r1smr, r0smr, lnrrsmr)) |>
    ungroup() -> fin


  ps1 <- hold %>%
    ungroup %>%
    dplyr::filter(sim == 1) %>%
    select(scenario, effect, pw, ps_results)

  list(ps_results = ps1,
       sim_ests = fin)
}

sum_results <- function(x){
  x %>%
  dplyr::group_by(scenario, effect, pw, name) |>
    dplyr::summarize(
      number_sims = n(),
      mean_value = mean(value),
      var_value = var(value)
    )
}

sum_rr <- function(x){
  x %>%
    dplyr::filter(grepl("rr", name)) %>%
    dplyr::mutate(
      name = gsub("ln", "", name),
      mean_value = exp(mean_value)
    )
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
      lnrr = log(r1/r0),
      n1ipt = sum(a*ipt),
      n0ipt = sum((1-a)*ipt),
      r1ipt = sum(a*y*ipt)/sum(a*ipt),
      r0ipt = sum((1-a)*y*ipt)/sum((1-a)*ipt),
      r1_mipt = (1/tot)*sum(a*y*ipt),
      r0_mipt = (1/tot)*sum((1-a)*y*ipt),
      lnrr_hajek_ipt = log(r1ipt/r0ipt),
      lnrr_ht_ipt = log(r1_mipt/r0_mipt),
      r1smr = sum(a*y*smr)/sum(a*smr),
      r0smr = sum((1-a)*y*smr)/sum((1-a)*smr),
      lnrrsmr = log(r1smr/r0smr)
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


