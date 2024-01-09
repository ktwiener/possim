library(dplyr)
bayes <- function(a, wa, w) {
  data.frame(
  pa1w1 = a*wa/w,
  pa1w0 = a*(1-wa)/(1-w)
  )
}

tibble::tribble(
  ~scenario, ~a, ~wa, ~w,
  "Complete", 0.2, 0.25, 0.25,
  "Complete", 0.2, 0.5, 0.5,
  "Complete", 0.2, 0.75, 0.75,
  "Partial",  0.2, 1, 0.25,
  "Partial",  0.2, 1, 0.5,
  "Partial",  0.2, 1, 0.75,
) %>%
  purrr::pmap_dfr(
    ~{
      dplyr::tibble(...) %>%
        dplyr::mutate(
          bayes(a, wa, w)
        )
    }
  )

outcome_prob <- function(a, w){
  expit(logit(0.01) + (logit(0.08)-logit(0.01))*w + a*log(0.8))
}


logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))

pa1w1 <- 0.2
pa1w0 <- 0.2

dplyr::tibble(
  scenario = rep(c("Complete", "Partial"), each = 3),
  pw = rep(c(0.25, 0.5, 0.75), times = 2),
  pa1w1 = c(0.2, 0.2, 0.2, 0.8, 0.4, 0.267),
  pa0w1 = 1-pa1w1,
  pa1w0 = c(0.2, 0.2, 0.2, 0, 0, 0),
  pa0w0 = 1-pa1w0
) %>%
  tidyr::pivot_longer(cols = c(pa1w1, pa0w1, pa1w0, pa0w0)) %>%
  dplyr::transmute(
    scenario, pw,
    a = if_else(grepl("a1", name), 1, 0),
    w = if_else(grepl("w1", name), 1, 0),
    conditional = value
  ) -> probs


dplyr::tibble(
  a = c(1, 1, 0, 0),
  w = c(1, 0, 1, 0),
  py = purrr::map2_dbl(a, w, outcome_prob)
) %>%
  dplyr::full_join(probs, by = c("a", "w")) %>%
  dplyr::mutate(
    pa_i = a*pa + (1-a)*(1-pa),
    pw_i = w*pw + (1-w)*(1-pw),
    paw  = conditional*pw_i
  ) %>%
  arrange(scenario, pw) %>%
  dplyr::mutate(pawy = paw*py) %>%
  dplyr::group_by(scenario, pw) %>%
  dplyr::summarize(py_marginal = sum(pawy))



