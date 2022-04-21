# BMGMR - Objective Bayesian meta-analysis 
This package implements the methods developed in [[1]](#1)

To reproduce the samples from the posterior using the methods developed in the paper you first run
```r
library(mvmeta)
set.seed(2021)

alp<-0.05
Np<-10^5 # number of observations drawn from the proposal distribution
B<-0.1*Np # burn in 
d<-3 ## degrees of freedom in t multivariate random effects model

# Determining X:p\times n and U:pn\times pn
dataREM<-hyp
# Matrix X
X<-t(cbind(dataREM$sbp,dataREM$dbp))
p<-nrow(X)  # model dimension
n<-ncol(X)  # sample size

# Matrix U
U<-matrix(0,n*p,n*p)
for (i_n in 1:n) {
    Use<-diag(c(dataREM$sbp_se[i_n],dataREM$dbp_se[i_n]))
    Corr_mat<-matrix(c(1,dataREM$rho[i_n],dataREM$rho[i_n],1),p,p)
    U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)]<- Use%*%Corr_mat%*%Use
}

# addtional definitons  
bi_n<-rep(1,n)
tbi_n<-t(bi_n)
In<-diag(bi_n)
Jn<-matrix(1,n,n)
Ip<-diag(rep(1,p))

# duplication matrix
Gp<-duplication_matrix(p)
Lp<-Gp%*%solve(t(Gp)%*%Gp)
```
To sample from the marginal posterior using a normal likelihood, using a jeffreys prior, you then do 
```r
res_nor_jefA<-sample_post_nor_jef_marg_mu(X,U,Np+B)
res_nor_jefB<-sample_post_nor_jef_marg_Psi(X,U,Np+B)
```
To sample from the marginal posterior using a normal likelihood, using the reference prior you use
```r
res_nor_refA<-sample_post_nor_ref_marg_mu(X,U,Np+B)
res_nor_refB<-sample_post_nor_ref_marg_Psi(X,U,Np+B)
```

The package also includes the `sample_post_t_jef_marg_mu`, `sample_post_t_jef_marg_Psi`, `sample_post_t_ref_marg_mu`and `sample_post_t_ref_marg_Psi` which does the same thing as the previous mentioned functions, though use the likelihood from the t-distribution. In order to make sure that the matrix `U` has the same covariance matrix, you need to standardize with the degrees of freedom. One example is
```r
Ut<-U*(d-2)/d
res_t_jefA<-sample_post_t_jef_marg_mu(X,Ut,d,Np+B)
```

## References
<a id="1">[1]</a> 
Olha Bodnar, Taras Bodnar (2021). 
Objective Bayesian meta-analysis based on generalized multivariate random effects model.
Under revision in Bayesian analysis.
