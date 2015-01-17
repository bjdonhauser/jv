# jv #

## Summary ##
* Package of R scripts for replication of results in Donhauser, B. J. (2012). "Jump Variation From High-Frequency Asset Returns: New Estimation Methods".
* Scripts related to the estimation of the empirical Bayesian jump variation extend and add to the functionality of the excellent [EbayesThresh](http://cran.r-project.org/web/packages/EbayesThresh/index.html) package.

## Usage Summary ##
* `bns04.R`: produces the jump variation estimator of Barndorff-Nielsen and Shephard (2004). This is the benchmark against which we compare our results. 
* `lm08.naive.R`: produces our naive estimator of jump variation based on jump location estimator of Lee and Mykland (2008). 
* `ebayesthresh.mod.R`: produces our empirical Bayesian estimators of jump variation.  
* `ret.cv.lm08.R`: provides a constant volatility simulation analysis of BNS04 and our naive estimator of jump variation.
* `ret.cv.R`: provides a constant volatility simulation analysis of BNS04 and our empirical Bayesian method of jump variation estimation.
* `ret.sv.R`: provides a stochastic volatility simulation analysis of BNS04 and our empirical Bayesian method of jump variation estimation.
* `ret.emp.bns04.R`: provides empirical analysis of BNS04 estimator of jump variation. 
* `ret.emp.lm08.R`: provides empirical analysis of our naive estimator of jump variation.
* `ret.emp.ebt.R`: provides empirical analysis of our empirical Bayesian method of jump variation estimation.
 

