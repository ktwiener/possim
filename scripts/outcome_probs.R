
library(dplyr)
library(furrr)
library(readr)
library(stringr)
library(rootSolve)

# Source all the R file with functions
all_scripts <- list.files("R", pattern = "*.R", full.names = T)
purrr::walk(all_scripts, source)

date <- "2024-05-02"
nsims <- 5000
positivity <- c("Full", "Partial")[1]
delt <- c("None","Homogeneous")[1]
probw <- c(25, 50, 75)[1]

scenarios <- expand.grid(positivity = positivity, delt = delt, probw = probw,
                         date = date, nsims = nsims, stringsAsFactors = FALSE)

outcome_prob <- function(positivity, delt, probw, date, nsims = 5000){
  patt <- paste0(tolower(c(positivity, delt, probw)), collapse = "-")
  files <- list.files(sprintf("data/population/%s", date), pattern = patt, full.names = T)
  pop <- read_csv(files)

  pop %>%
    summarize(
      scenario=scenario[1],
      effect = effect[1],
      pw = pw[1],
      delta = delta[1],
      outcome_prop = sum(y)/n()
    )
}

marginal_outcome_probabilities <- purrr::pmap_dfr(scenarios, outcome_prob)
