
simulation_box <- function(x, eff = c("None", "Homogeneous"), hajek = T, prefix = "mest"){
  keep <- grepl("hajek", x$pars)
  if (!hajek) {
    keep <- !keep
    filenm <- "ht"
  } else filenm <- "hajek"

  keep <- keep | grepl("smr", x$pars)

  truth <- x$delta[x$effect == eff][1]
    fig <-
    x %>%
    select(scenario, effect, pw, sims, delta, pars, ests) %>%
    dplyr::mutate(scenario = gsub(" exchangeability", "", scenario),
                  scenario = gsub("Full", "Complete", scenario)) %>%
    filter(effect == eff, keep) %>%
    dplyr::mutate(pars = factor(pars,
                                levels = c("ef_ipt_hajek_lnrr", "lnrr_hajek_ipt", "ef_ipt_lnrr", "lnrr_ht_ipt", "ef_smr_lnrr", "lnrrsmr"),
                                labels = c("IPT (Hajek)","IPT (Hajek)", "IPT (Horvitz-Thompson)", "IPT (Horvitz-Thompson)", "SMR", "SMR"))) %>%
    dplyr::mutate(pw = paste0("P(W=1) = ",pw),
                  bias = ests - delta) %>%
    ggplot(mapping = aes(y = bias, x = pars, color = pars)) +
    geom_hline(yintercept = 0) +
    geom_boxplot(coef = 0, outlier.shape = NA, varwidth = T) +
    scale_y_continuous(limits = c(-1, 2)) +
    geom_violin(alpha = 0.5, aes(fill = pars)) +
    facet_grid(vars(scenario), vars(pw), switch = "y") +
    scale_color_grey(guide = "none") +
    scale_fill_grey() +
    theme_classic(base_size = 12) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(margin = margin(r = -50), size = 12),
          strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"),
          strip.text.x.top = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.spacing.x = unit(0, 'lines'),
          plot.margin = unit(c(0, 0, 0, 0.6),
                             "inches"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = "bold")
    ) +
    ylab("Error") +
    labs(fill = "Weight")

  ggsave(filename = sprintf("./data/results/figures/%s-%s-%s-plot.png", prefix, tolower(eff), filenm),
         plot = fig,
         units = "in",
         width = 5.5,
         height = 6.5
         )

  return(fig)

}

simulation_violin <- function(x, eff = c("None", "Homogeneous"), hajek = T, prefix = "mest"){
  keep <- grepl("hajek", x$pars)
  if (!hajek) {
    keep <- !keep
    filenm <- "ht"
  } else filenm <- "hajek"

  keep <- keep | grepl("smr", x$pars)

  truth <- x$delta[x$effect == eff][1]
  fig <-  x %>%
    select(scenario, effect, pw, sims, delta, pars, ests) %>%
    dplyr::mutate(scenario = gsub(" exchangeability", "", scenario)) %>%
    filter(effect == eff, keep) %>%
    dplyr::mutate(pars = factor(pars,
                                levels = c("ef_ipt_hajek_lnrr", "ef_ipt_lnrr", "ef_smr_lnrr"),
                                labels = c("IPT (Hajek)", "IPT (Horvitz-Thompson)", "SMR"))) %>%
    dplyr::mutate(pw = paste0("P(W=1) = ",pw)) %>%
    ggplot(mapping = aes(y = ests, x = pars, color = pars, fill = pars)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = truth, alpha = 0.5, linetype = 2) +
    geom_violin(alpha = 0.5) +
    scale_y_continuous(limits = c(-1.6, 1.6)) +
    facet_grid(vars(scenario), vars(pw), switch = "y") +
    scale_color_manual(values=wes_palette(n=2, name="Royal1"), guide = "none") +
    scale_fill_manual(values=wes_palette(n=2, name="Royal1")) +
    theme_classic(base_size = 12) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(margin = margin(r = -55), size = 11),
          strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"),
          strip.text.x.top = element_text(size = 12, face = "bold"),
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.spacing.x = unit(0, 'lines'),
          plot.margin = unit(c(0, 0, 0, 0.6),
                             "inches"),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12, face = "bold")
    ) +
    ylab("log RR") +
    labs(fill = "Estimator")

  ggsave(filename = sprintf("./data/results/%s-%s-%s-violin.png", prefix, tolower(eff), filenm),
         plot = fig,
         units = "in",
         width = 5.5,
         height = 6.5)

  return(fig)

}
