#' Interface for the Metropolis-Hasting algorithms
#'
#' @param X
#' @param U
#' @param Np number of samples to simulate from the markov chain
#' @param burn_in the number of burnin samples
#'
#' @export
BMGMR <- function(X, U, N, burn_in, likelihood, prior, algorithm_version, d=NULL) {
  assertthat::assert_that(likelihood %in% c("normal", "t"))
  assertthat::assert_that(prior %in% c("jeffrey", "reference"))
  assertthat::assert_that(algorithm_version %in% c("A", "B"))

  if (likelihood == "normal" & prior == "jeffrey" & algorithm_version == "A") {
    simulations <- sample_post_nor_jef_marg_mu(X, U, N + burn_in)
  }else if (likelihood == "normal" & prior == "reference" & algorithm_version == "A") {
    simulations <- sample_post_nor_ref_marg_mu(X, U, N + burn_in)
  }else if (likelihood == "t" & prior == "jeffrey") {
    simulations <-
  }else if (likelihood == "t" & prior == "jeffrey") {
  }
}

#' Summary statistics from a posterior distribution
#'
#' @param x a sample from the posterior which to aggregate from.
#' @param alp the confidence level used in computation of the summary statistics
#'
#' @return a vector with summary statistics
#' @export
Bayes_inference <- function(x,alp){
  x_mean<-apply(x,1,mean)
  x_med<-apply(x,1,quantile,probs=0.5)
  x_ql<-apply(x,1,quantile,probs=alp/2)
  x_qu<-apply(x,1,quantile,probs=1-alp/2)
  x_sd<-sqrt(apply(x,1,var))
  rbind(x_mean,x_med,x_sd,x_ql,x_qu)
}
