# BMGMR - Objective Bayesian meta-analysis 
This package implements the methods developed in [[1]](#1)

To simulate samples from the posterior using the methods developed in the paper you simply need to run
```r

set.seed(2021)
dataREM<-mvmeta::hyp
# Observation matrix X
X<-t(cbind(dataREM$sbp,dataREM$dbp))
p<-nrow(X) # model dimension
n<-ncol(X) # sample size
# Matrix U
U<-matrix(0,n*p,n*p)
for (i_n in 1:n) {
  Use<-diag(c(dataREM$sbp_se[i_n],dataREM$dbp_se[i_n]))
  Corr_mat<-matrix(c(1,dataREM$rho[i_n],dataREM$rho[i_n],1),p,p)
  U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)]<- Use%*%Corr_mat%*%Use
}
bmgmr_run <- BMGMR(X, U, 1e4, burn_in = 100,
                   likelihood = "normal", prior="jeffrey",
                   algorithm_version = "A")
summary(bmgmr_run)
plot(bmgmr_run, range.x=list(c(-14,-6), c(-7,-3)), nbins = 51, 
     ci.levels=c(0.90,0.95,0.99),col=topo.colors(5)[-1])
```

## References
<a id="1">[1]</a> 
Olha Bodnar, Taras Bodnar (2021). 
Objective Bayesian meta-analysis based on generalized multivariate random effects model.
Under revision in Bayesian analysis.
