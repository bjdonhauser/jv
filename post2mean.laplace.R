"post2mean.laplace" <-
function(x, w, a = 0.5)
{
#
# find the posterior squared mean for the double exponential prior for 
#   given x, w and a, assuming the error variance
#   is 1.
#
  #  only allows a < 20.
  a <- min(a, 20)

  #  First find the odds of zero and the shrinkage factor
	wpost <- 1 - (1 - w)/(1 + w * beta.laplace(x, a))	
  
	#  now find the posterior mean conditional on being non-zero
	sx   <- sign(x)
	x    <- abs(x)
  xpa  <- x + a
  xma  <- x - a
	cp1  <- pnorm(xma)
	dp1  <- dnorm(xma)
	cp2  <- pnorm(-xpa)
	dp2  <- dnorm(xpa)
  den  <- exp(-a*x)*cp1 + exp(a*x)*cp2
  cnum <- exp(-a*x)*cp1 - exp(a*x)*cp2
  dnum <- exp(-a*x)*dp1 + exp(a*x)*dp2
	#ef <- exp(pmin(2 * a * x, 100))  
  # legacy from postmean.laplace.R
  # does implicit robustification of x
  
  # Including an approximation. 
  # Unstable
#   postmeancond <- xma*xma + 1  
#   postmeancond[xpa < 6] <- x[xpa < 6]^2 + a^2 + 1 - 
#                            2*a*x[xpa < 6] * cnum[xpa < 6] / den[xpa < 6] - 
#                            2*a * exp(-a*x[xpa < 6]) dp1[xpa < 6] / 
#                            den[xpa < 6]
  postmeancond <- x^2 + a^2 + 1 - 2*a*x*cnum/den - 2*a*exp(-a*x)*dp1/den 
                               
#  calculate posterior mean and return
#
	mu2tilde <- wpost * postmeancond  # no `sx` as in `postmean.R`
	return(mu2tilde)
}