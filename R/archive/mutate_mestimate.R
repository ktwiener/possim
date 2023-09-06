ests <- pop |>
  group_by(scenario, delta, sim) |>
  tidyr::nest(data = c(3, 5:last_col())) |>
  dplyr::mutate(
    ests = purrr::map(data,
                      function(dset) {
                        m <- geex::m_estimate(
                          estFUN = eefun_ipt,
                          data = as.data.frame(dset),
                          root_control = setup_root_control(FUN = lm_redo,
                                                            start = c(-2, .5, rep(.5, 8)),
                                                            roots_name = "par")
                        )

                        list(
                          pars = c("ps_beta0", "ps_beta1",
                                   "ef_ipt_r1", "ef_ipt_r0", "ef_ipt_lnrr", "ef_ipt_lnor",
                                   "ef_smr_r1", "ef_smr_r0", "ef_smr_lnrr", "ef_smr_lnor"),
                          ests = m@estimates,
                          vars = diag(m@vcov)
                        )
                      }
    )
  )
