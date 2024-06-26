
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
positivity <- c("Full", "Partial")
delt <- c("None","Homogeneous")
probw <- c(25, 50, 75)

scenarios <- expand.grid(positivity = positivity, delt = delt, probw = probw,
                         date = date, nsims = nsims, stringsAsFactors = FALSE)

run_sims <- function(positivity, delt, probw, date, nsims = 5000){
  patt <- paste0(tolower(c(positivity, delt, probw)), collapse = "-")
  files <- list.files(sprintf("data/population/%s", date), pattern = patt, full.names = T)
  pop <- read_csv(files)

  h <- pop %>%
    dplyr::group_by(sim) %>%
    tidyr::nest()

  plan(multisession(workers = 7))
  res <- furrr::future_map_dfr(h$data[1:nsims], mestimate)

  saveRDS(res, sprintf("data/simulations/%s/%s-%s.rds", date, patt, nsims))
}

s <- Sys.time()
purrr::pwalk(
  scenarios,
  run_sims
)
fin_fut <- Sys.time() - s
