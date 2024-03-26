each_pop <- function(popnm, done_files, sim_loc, sim_number = 1:5000){

  type <- stringr::str_split(popnm, "/|\\.")[[1]][4]

  done_sims <- as.numeric(
    substr(done_files[grepl(type, done_files)], nchar(done_files)-8, nchar(done_files)-4))

  pop <- readRDS(popnm)

  prep_pop <- pop |>
    dplyr::filter(sim %in% sim_number, !sim %in% done_sims) %>%
    group_by(scenario, effect, pw, sim) |>
    tidyr::nest() %>%
    dplyr::mutate(
      filename = tolower(sprintf("%s-%s-%s-%05g", word(scenario, 1), effect,
                                 as.numeric(pw)*100, sim))
    )


  if (nrow(prep_pop) == 0) return()

  future::plan(multisession)
  start <- Sys.time()
  furrr::future_walk2(prep_pop$data,
                      prep_pop$filename,
                      function(dset, filenm) {
                        m <- geex::m_estimate(
                          estFUN = eefun_ipt,
                          data = as.data.frame(dset),
                          root_control = setup_root_control(FUN = lm_redo,
                                                            start = c(-2, .5, rep(.5, 4*3)),
                                                            roots_name = "par"))
                        ret <- list(
                          pars = c("ps_beta0", "ps_beta1",
                                   "ef_ipt_r1", "ef_ipt_r0", "ef_ipt_lnrr", "ef_ipt_lnor",
                                   "ef_ipt_hajek_r1", "ef_ipt_hajek_r0", "ef_ipt_hajek_lnrr", "ef_ipt_hajek_lnor",
                                   "ef_smr_r1", "ef_smr_r0", "ef_smr_lnrr", "ef_smr_lnor"),
                          ests = m@estimates,
                          vars = diag(m@vcov)
                        )

                        saveRDS(ret, sprintf("%s/%s.rds", sim_loc, filenm))
                      }
  )

  end <- Sys.time()
  furrr_time2 <- end - start
}
