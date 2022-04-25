#' Interface for the BMGMR class
#'
#' This is the main interface to the MCMC algorithms implemented in this
#' package. The data is on wide format with p rows (variables) and n columns
#' (observations).
#'
#' @param X observation matrix of dimension p x n
#' @param U block matrix with uncertainties of size pn x pn
#' @param N number of samples to simulate from the markov chain
#' @param burn_in the number of burn-in samples
#' @param likelihood the likelihood to use. It currently supports "normal" and
#' "t".
#' @param prior the prior to use. It currently supports "reference" and
#' "jeffrey".
#' @param algorithm_version One of "A" or "B". Both algorithms samples the same
#' quantities.
#' @param d the degrees of freedom which is used when the "t" option is used for
#' the likelihood.
#'
#' @return a BMGMR class.
#' @export
BMGMR <- function(X, U, N, burn_in, likelihood, prior, algorithm_version, d=NULL) {
  assertthat::assert_that(likelihood %in% c("normal", "t"))
  assertthat::assert_that(prior %in% c("jeffrey", "reference"))
  assertthat::assert_that(algorithm_version %in% c("A", "B"))

  if (algorithm_version == "A") {
    if (likelihood == "normal" & prior == "jeffrey") {
      simulations <- sample_post_nor_jef_marg_mu(X, U, N + burn_in)
    }else if (likelihood == "normal" & prior == "reference") {
      simulations <- sample_post_nor_ref_marg_mu(X, U, N + burn_in)
    }else if (likelihood == "t" & prior == "jeffrey") {
      simulations <- sample_post_t_jef_marg_mu(X, U*(d-2)/d, d, N + burn_in)
    }else if (likelihood == "t" & prior == "reference") {
      simulations <- sample_post_t_ref_marg_mu(X, U*(d-2)/d, d, N + burn_in)
    }
  }else{
    if (likelihood == "normal" & prior == "jeffrey") {
      simulations <- sample_post_nor_jef_marg_Psi(X, U, N + burn_in)
    }else if (likelihood == "normal" & prior == "reference") {
      simulations <- sample_post_nor_ref_marg_Psi(X, U, N + burn_in)
    }else if (likelihood == "t" & prior == "jeffrey") {
      simulations <- sample_post_t_jef_marg_Psi(X, U*(d-2)/d, d, N + burn_in)
    }else if (likelihood == "t" & prior == "reference") {
      simulations <- sample_post_t_ref_marg_Psi(X, U*(d-2)/d, d, N + burn_in)
    }
  }

  structure(list(mu=simulations[[1]],
       psi=simulations[[2]],
       X=X,
       U=U,
       N=N,
       p=nrow(X),
       burn_in=burn_in,
       likelihood=likelihood,
       prior=prior,
       algorithm_version=algorithm_version,
       d=d
       ), class="BMGMR")
}

#' Summary statistics from the posterior of a BMGMR class
#'
#' It will use the data without burn in
#'
#' @param object is a BMGMR class
#' @param alpha is the confidence level used in construction of the credible
#' intervals
#' @param ... not used
#'
#' @returns a list with summary statistics
#' @export
summary.BMGMR <- function(object, alpha=0.95, ...) {
  Gp<-duplication_matrix(object$p)
  Lp<-Gp%*%solve(t(Gp)%*%Gp)

  list("mu"=bayes_inference(object$mu[,(object$burn_in+1):(object$burn_in+object$N)], alpha),
       "psi"=bayes_inference(object$psi[,(object$burn_in+1):(object$burn_in+object$N)], alpha)%*%Lp)
}

#' Summary statistics from a posterior distribution
#'
#' @param x a sample from the posterior which to aggregate from.
#' @param alp the confidence level used in computation of the summary statistics
#'
#' @return a matrix with summary statistics
bayes_inference <- function(x,alp){
  x_mean<-apply(x,1,mean)
  x_med<-apply(x,1,quantile,probs=0.5)
  x_ql<-apply(x,1,quantile,probs=alp/2)
  x_qu<-apply(x,1,quantile,probs=1-alp/2)
  x_sd<-sqrt(apply(x,1,var))
  rbind(x_mean,x_med,x_sd,x_ql,x_qu)
}


#' Plot a BMGMR object
#'
#' @param x a BMGMR object
#' @param ... optional arguments to ci2d (if gplots is installed)
#' @export
plot.BMGMR <- function(x, ...) {
  for (var in c("mu", "psi")) {
    df <- x[[var]]
    ylab <- function(p) ifelse(var == "mu", parse(text=paste0("mu[", p, "]")),
                               parse(text=paste0("psi[", p, "]")))
    for (idx in 1:nrow(df)) {
      plot(y=df[idx,], x=1:ncol(df), type="l", ylab=ylab(idx),
           xlab="index")
      abline(v=x$burn_in, col = "lightgray", lty = 3)
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }

  if (!requireNamespace("gplots")) {
    print("Could not load the gplots package.")
  }else{
    gplots::ci2d(t(x[["mu"]]),...)
    title(main=paste0("Model:", x$likelihood,"-", x$prior, ", Algorithm ",
                      x$algorithm_version), xlab="SBP",ylab="DBP")
  }
}
