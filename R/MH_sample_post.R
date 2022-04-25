
#' A Duplication matrix
#'
#' This function creates a duplication matrix of size x
#'
#' @param x the dimension of the matrix
#'
#' @return a matrix of size x^2
duplication_matrix <- function(x){
  mat <- diag(x)
  index <- seq(x*(x+1)/2)
  mat[ lower.tri( mat , TRUE ) ] <- index
  mat[ upper.tri( mat ) ] <- t( mat )[ upper.tri( mat ) ]
  outer(c(mat), index , function( x , y ) ifelse(x==y, 1, 0 ) )
}

#' Metropolis hasting algorithm A
#'
#' This function implements the algorithm A using the normal likelihood and
#' jeffreys prior. The number of observations is n and the number of variables
#' is p.
#'
#' @param X A p x n matrix which contains the observations
#' @param U A p n x p n block matrix which contains uncertainties.
#' @param Np the number of simulations to perform.
#'
#' @return list with samples from the two marginal distributions
sample_post_nor_jef_marg_mu<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size

  ############## addtional definitons
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-duplication_matrix(p)
  tGp<-t(Gp)

  mu_m<-NULL
  Psi_m<-NULL

  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))

  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n+p),p,n+1)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p

  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+2))*exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))

  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
  }
  p0<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num

  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal
    mu_p<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n+p),p,n+1)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p

    ### value of the proposal at new draw
    q1 <-det(Psi_p)^(-0.5*(n+p+2))*exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))

    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    }
    p1<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num

    ### MH ratio
    ratio_MH<-p1*q0/p0/q1

  if (runif(1) <= ratio_MH)
    {
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }

    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
  }
  output<-list(mu_m,Psi_m)
}

#' Metropolis Hasting Algorithm B, normal-distribution
#'
#' A MH algorithm that samples using a likelihood from the normal distribution
#' and the reference prior.
#'
#' @inherit sample_post_nor_jef_marg_mu
sample_post_nor_jef_marg_Psi<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size

  ############## addtional definitons
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-duplication_matrix(p)
  tGp<-t(Gp)

  mu_m<-NULL
  Psi_m<-NULL

  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  cS<-sqrt(n-1)*chol(S)

  ### generating an initial value for mu and Psi new draw from proposal
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cS)%*%solve(Z%*%t(Z))%*%cS
  mu0<-bar_X+t(chol(Psi0))%*%rnorm(p)/sqrt(n)
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)

  ### initial value of the proposal
  exp_prop<-exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
  q0<-det(Psi0)^(-0.5*(n+p+2))*exp_prop

  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  #exp_prop<-1
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
  }
  p0<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num

  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cS)%*%solve(Z%*%t(Z))%*%cS
    mu_p<-bar_X+t(chol(Psi_p))%*%rnorm(p)/sqrt(n)
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)

    ### value of the proposal at new draw
    exp_prop<-exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))
    q1<-det(Psi_p)^(-0.5*(n+p+2))*exp_prop

    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    #exp_prop<-1
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    }
    p1<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num

    ### MH ratio
    ratio_MH<-p1*q0/p0/q1

    if (runif(1) <= ratio_MH)
    {
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }

    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
  }
  output<-list(mu_m,Psi_m)
}

#' Metropolis Hasting Algorithm A, t-distribution
#'
#' A MH algorithm that samples using a likelihood from the t distribution and
#' the jeffreys prior.
#'
#' @inherit sample_post_nor_jef_marg_mu
#' @param d the degrees of freedom for the t-distribution
#'
sample_post_t_jef_marg_mu<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size

  ############## addtional definitons
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-duplication_matrix(p)
  tGp<-t(Gp)

  mu_m<-NULL
  Psi_m<-NULL

  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))

  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n+p),p,n+1)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))

  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))

  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  t_num<-0
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
  }
  p0<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))

  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal
    mu_p<-bar_X+sqrt((n-1)/n/(n-p+1))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p+1)/(n-p+1))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n+p),p,n+1)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))

    ### value of the proposal at new draw
    q1 <-det(Psi_p)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))

    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    t_num<-0
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    }
    p1<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))

    ### MH ratio
    ratio_MH<-p1*q0/p0/q1

    if (runif(1) <= ratio_MH)
    {
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }

    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
  }
  output<-list(mu_m,Psi_m)
}

#' Metropolis Hasting Algorithm B, t-distribution
#'
#' A MH algorithm that samples using a likelihood from the t distribution
#' and the jeffreys prior.
#'
#' @inherit sample_post_t_jef_marg_mu
sample_post_t_jef_marg_Psi<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size

  ############## addtional definitons
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-duplication_matrix(p)
  tGp<-t(Gp)

  mu_m<-NULL
  Psi_m<-NULL

  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  cS<-sqrt(n-1)*chol(S)

  ### generating an initial value for mu and Psi new draw from proposal
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cS)%*%solve(Z%*%t(Z))%*%cS*sqrt(rchisq(1,d)/(d))
  mu0<-bar_X+sqrt((d+(n-1)*sum(diag(solve(Psi0)%*%S)))/n/(d+p*n-p))*t(chol(Psi0))%*%rnorm(p)/sqrt(rchisq(1,d+p*n-p)/(d+p*n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)

  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))

  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  t_num<-0
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
  }
  p0<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))

  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cS)%*%solve(Z%*%t(Z))%*%cS*sqrt(rchisq(1,d)/(d))
    mu_p<-bar_X+sqrt((d+(n-1)*sum(diag(solve(Psi_p)%*%S)))/n/(d+p*n-p))*t(chol(Psi_p))%*%rnorm(p)/sqrt(rchisq(1,d+p*n-p)/(d+p*n-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)

    ### value of the proposal at new draw
    q1<-det(Psi_p)^(-0.5*(n+p+2))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))

    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    t_num<-0
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    }
    p1<-sqrt(det(mat_numJ/n))*sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))

    ### MH ratio
    ratio_MH<-p1*q0/p0/q1

    if (runif(1) <= ratio_MH)
    {
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }

    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
  }
  output<-list(mu_m,Psi_m)
}


#' Metropolis Hasting Algorithm A, normal-distribution
#'
#' A MH algorithm that samples using a likelihood from the normal distribution
#' and the reference prior.
#'
#' @inherit sample_post_nor_jef_marg_mu
sample_post_nor_ref_marg_mu<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size

  ############## addtional definitons
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-duplication_matrix(p)
  tGp<-t(Gp)

  mu_m<-NULL
  Psi_m<-NULL

  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))

  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p

  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+1))*exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))

  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  det_num<-det_num*sqrt(det(iPsiUi))
  exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
  }
  p0<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num

  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal
    mu_p<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p

    ### value of the proposal at new draw
    q1 <-det(Psi_p)^(-0.5*(n+p+1))*exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))

    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    }
    p1<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num

    ### MH ratio
    ratio_MH<-p1*q0/p0/q1

    if (runif(1) <= ratio_MH)
    {mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }

    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
  }
  output<-list(mu_m,Psi_m)
}


#' Metropolis Hasting Algorithm B, normal-reference
#'
#' A MH algorithm that samples using a likelihood from the normal distribution
#' and the reference prior.
#'
#' @inherit sample_post_nor_jef_marg_mu
sample_post_nor_ref_marg_Psi<-function(X,U,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size

  ############## addtional definitons
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-duplication_matrix(p)
  tGp<-t(Gp)

  mu_m<-NULL
  Psi_m<-NULL

  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  cS<-sqrt(n-1)*chol(S)

  ### generating an initial value for mu and Psi new draw from proposal
  Z<-matrix(rnorm(p*n-p),p,n-1)
  Psi0<-t(cS)%*%solve(Z%*%t(Z))%*%cS
  mu0<-bar_X+t(chol(Psi0))%*%rnorm(p)/sqrt(n)
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)

  ### initial value of the proposal
  exp_prop<-exp(-0.5*sum(diag(solve(Psi0)%*%Cov_p)))
  q0<-det(Psi0)^(-0.5*(n+p+1))*exp_prop

  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  #exp_prop<-1
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  det_num<-det_num*sqrt(det(iPsiUi))
  exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
  }
  p0<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num

  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal
    Z<-matrix(rnorm(p*n-p),p,n-1)
    Psi_p<-t(cS)%*%solve(Z%*%t(Z))%*%cS
    mu_p<-bar_X+t(chol(Psi_p))%*%rnorm(p)/sqrt(n)
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)

    ### value of the proposal at new draw
    exp_prop<-exp(-0.5*sum(diag(solve(Psi_p)%*%Cov_p)))
    q1<-det(Psi_p)^(-0.5*(n+p+1))*exp_prop

    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    #exp_prop<-1
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    det_num<-det_num*sqrt(det(iPsiUi))
    exp_num<-exp_num*exp(-0.5*sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n]))
    }
    p1<-sqrt(det(tGp%*%mat_num%*%Gp/n))*det_num*exp_num

    ### MH ratio
    ratio_MH<-p1*q0/p0/q1

    if (runif(1) <= ratio_MH)
    {
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }

    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
  }
  output<-list(mu_m,Psi_m)
}

#' Metropolis Hasting Algorithm A, t-reference
#'
#' A MH algorithm that samples using a likelihood from the t distribution
#' and the reference prior.
#'
#' @inherit sample_post_t_jef_marg_mu
sample_post_t_ref_marg_mu<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size

  ############## addtional definitons
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-duplication_matrix(p)
  tGp<-t(Gp)

  mu_m<-NULL
  Psi_m<-NULL

  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  tcholS<-t(chol(S))

  ### generating an initial value for mu and Psi new draw from proposal
  mu0<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)
  cCov_p<-chol(Cov_p)
  Z<-matrix(rnorm(p*n),p,n)
  Psi0<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))

  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))

  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  t_num<-0
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
  }
  p0<-sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))

  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal
    mu_p<-bar_X+sqrt((n-1)/n/(n-p))*tcholS%*%rnorm(p)/sqrt(rchisq(1,n-p)/(n-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)
    cCov_p<-chol(Cov_p)
    Z<-matrix(rnorm(p*n),p,n)
    Psi_p<-t(cCov_p)%*%solve(Z%*%t(Z))%*%cCov_p*sqrt(rchisq(1,d)/(d))

    ### value of the proposal at new draw
    q1 <-det(Psi_p)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))

    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    t_num<-0
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    }
    p1<-sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))

    ### MH ratio
    ratio_MH<-p1*q0/p0/q1

    if (runif(1) <= ratio_MH)
    {
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }

    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
  }
  output<-list(mu_m,Psi_m)
}

#' Metropolis Hasting Algorithm B, t-reference
#'
#' A MH algorithm that samples using a likelihood from the t distribution
#' and the reference prior.
#'
#' @inherit sample_post_t_jef_marg_mu
sample_post_t_ref_marg_Psi<-function(X,U,d,Np){
  p<-nrow(X)  # model dimension
  n<-ncol(X)  # sample size

  ############## addtional definitons
  bi_n<-rep(1,n)
  tbi_n<-t(bi_n)
  In<-diag(bi_n)
  Jn<-matrix(1,n,n)
  Ip<-diag(rep(1,p))
  Gp<-duplication_matrix(p)
  tGp<-t(Gp)

  mu_m<-NULL
  Psi_m<-NULL

  bar_X<-X%*%bi_n/n
  S<-X%*%(In-Jn/n)%*%t(X)/(n-1)
  cS<-sqrt(n-1)*chol(S)

  ### generating an initial value for mu and Psi new draw from proposal
  Z<-matrix(rnorm(p*(n-1)),p,n-1)
  Psi0<-t(cS)%*%solve(Z%*%t(Z))%*%cS*sqrt(rchisq(1,d)/(d))
  mu0<-bar_X+sqrt((d+(n-1)*sum(diag(solve(Psi0)%*%S)))/n/(d+p*n-p))*t(chol(Psi0))%*%rnorm(p)/sqrt(rchisq(1,d+p*n-p)/(d+p*n-p))
  Xn<-X-mu0%*%tbi_n
  Cov_p<-Xn%*%t(Xn)

  ### initial value of the proposal
  q0<-det(Psi0)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi0)%*%Cov_p))/d)^(-0.5*(p*n+d))

  ### initial value of the target (posterior)
  mat_num<-matrix(0,p^2,p^2)
  mat_numJ<-matrix(0,p,p)
  det_num<-1
  exp_num<-1
  t_num<-0
  for (i_n in 1:n)
  {iPsiUi<-solve(Psi0+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
  mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
  mat_numJ<-mat_numJ+iPsiUi
  det_num<-det_num*sqrt(det(iPsiUi))
  t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
  }
  p0<-sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))

  for (j in 1:Np)
  {
    ### generating a new draw mu and Psi from the proposal
    Z<-matrix(rnorm(p*(n-1)),p,n-1)
    Psi_p<-t(cS)%*%solve(Z%*%t(Z))%*%cS*sqrt(rchisq(1,d)/(d))
    mu_p<-bar_X+sqrt((d+(n-1)*sum(diag(solve(Psi_p)%*%S)))/n/(d+p*n-p))*t(chol(Psi_p))%*%rnorm(p)/sqrt(rchisq(1,d+p*n-p)/(d+p*n-p))
    Xn<-X-mu_p%*%tbi_n
    Cov_p<-Xn%*%t(Xn)

    ### value of the proposal at new draw
    q1<-det(Psi_p)^(-0.5*(n+p+1))*(1+sum(diag(solve(Psi_p)%*%Cov_p))/d)^(-0.5*(p*n+d))

    ### value of the posterior
    mat_num<-matrix(0,p^2,p^2)
    mat_numJ<-matrix(0,p,p)
    det_num<-1
    exp_num<-1
    t_num<-0
    for (i_n in 1:n)
    {iPsiUi<-solve(Psi_p+U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)])
    mat_num<-mat_num+kronecker(iPsiUi,iPsiUi)
    mat_numJ<-mat_numJ+iPsiUi
    det_num<-det_num*sqrt(det(iPsiUi))
    t_num<-t_num+sum(t(Xn[,i_n])%*%iPsiUi%*%Xn[,i_n])/d
    }
    p1<-sqrt(det(tGp%*%(mat_num-as.vector(mat_numJ)%*%t(as.vector(mat_numJ))/(n*p+d))%*%Gp/n))*det_num*(1+t_num)^(-0.5*(p*n+d))

    ### MH ratio
    ratio_MH<-p1*q0/p0/q1

    if (runif(1) <= ratio_MH)
    {
    mu0<-mu_p
    Psi0<-Psi_p
    q0<-q1
    p0<-p1
    }

    Psi_m<-cbind(Psi_m,as.vector(Psi0))
    mu_m<-cbind(mu_m,mu0)
  }
  output<-list(mu_m,Psi_m)
}
