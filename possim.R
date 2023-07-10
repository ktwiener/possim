## Script to run through simulations

devtools::load_all()

pop <- purrr::pmap_dfr(
  .l = settings,
  .f = popgen
)

pspop <- pop |>
  group_by(sim) |>
  psest()

dplyr::group_modify(.data = pspop,
                 .f = ~cond_effects(.x))

dplyr::group_modify(.data = pspop,
                    .f = ~ true_effects(.x))
