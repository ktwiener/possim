# Helper functions and utilities
#
# This file contains small functions that come in handy in several
# places, or are a bit too general to put in other scripts.
#

#' Expit function
#'
#' @description Takes the inverse of the log odds to return a probability.
#'
#' @param x A number (preferably a log-odds)
#' @return A probability.
#' @export


expit <- function(x) exp(x)/(1+exp(x))

logit <- function(p)log(p/(1-p))

or_func <- function(risk1, risk0) risk1*(1-risk0)/(risk0*(1-risk1))
#' The Bernoulli distribution
#'
#' @description Random generation for the Bernoulli distribution with parameters `prob`
#'
#' @param n number of observations
#' @param prob probability of success on each trial
#' @return vector of length n
#' @export

rbern <- function(n, prob) rbinom(n, 1, prob)

## Read in individual datasets
read_chunked <- function(filepat){
  all_sims <- list.files(path = "./data/simulations/", pattern = filepat, full.names = T)
  purrr::map_dfr(
    all_sims,
    readRDS
  )
}
