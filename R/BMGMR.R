


#' Interface for the Metropolis-Hasting algorithms
#'
#'
#' @export
BMGMR <- function(X,U,Np,d=NULL) {

}

#' Summary statistics from a posterior distribution
#'
#' @param x a sample from the posterior which to aggregate from.
Bayes_inference <- function(x,alp){
  x_mean<-apply(x,1,mean)
  x_med<-apply(x,1,quantile,probs=0.5)
  x_ql<-apply(x,1,quantile,probs=alp/2)
  x_qu<-apply(x,1,quantile,probs=1-alp/2)
  x_sd<-sqrt(apply(x,1,var))
  rbind(x_mean,x_med,x_sd,x_ql,x_qu)
}
