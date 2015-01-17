"ebayesthresh.mod" <-
function(x, prior = "laplace", a = 0.5, bayesfac = FALSE, sdev = NA, verbose = FALSE, 
	threshrule = "median", window.type = "full", window.len = NA)
{
#  Given a vector of data x, find the marginal maximum likelihood estimator
#   of the mixing weight w, and apply an appropriate thresholding rule using
#   this weight.
#  If the prior is laplace and a=NA, then the scale factor is also found by MML.
#  The thresholding rules allowed are "median", "mean", "hard", "soft" and "none";
#   if "none" is used, then only the parameters are worked out.
#  If hard or soft thresholding is used, the argument "bayesfac" specifies
#   whether to use the bayes factor threshold or the posterior median threshold.
#  If verbose=TRUE then the routine returns a list with several arguments, including
#   muhat which is the result of the thresholding.
#  If verbose=FALSE then only muhat is returned.
#  It is assumed that the standard deviation of the data is sdev; if sdev=NA, then
#   it is estimated using the function mad(x).
#
#  find the standard deviation if necessary and estimate the parameters
#
#  Requires:
#    wfromx.mod(), tfromw.mod(), post2mean.laplace(), jv.median()
#
	sdev <- lm08.sdev(x, sdev, window.type, window.len)
  #cat("sdev is: ", sdev)
  sdev[sdev == 0] <- min(sdev[sdev != 0])  # a concern when window.len=5
	#cat("sdev still is: ", sdev)
	x <- x/sdev
	pr <- substring(prior, 1, 1)
	if((pr == "l") & is.na(a)) {
		pp <- wandafromx(x)
		w <- pp$w
		a <- pp$a
	}
	else w <- wfromx.mod(x, prior = prior, a = a)	#
	if(pr != "m" | verbose) {
		tt <- tfromw.mod(w, prior = prior, bayesfac = bayesfac, a = a)
		tcor <- sdev * tt
	}
	if(threshrule == "median")
		muhat <- postmed(x, w, prior = prior, a = a)
	if(threshrule == "mean")
		muhat <- postmean(x, w, prior = prior, a = a)
	if(threshrule == "hard")
		muhat <- threshld(x, tt)
	if(threshrule == "soft")
		muhat <- threshld(x, tt, hard = FALSE)
	if(threshrule == "none") muhat <- NA	#
   
  # Now build desired output 
  muhat      <- sdev * muhat
  muhat2     <- muhat^2
  jvhat.js05 <- sum(muhat2)
 
  # my additions: mu2hat and jvhat
  if (pr == "l") {
    mu2hat <- post2mean.laplace(x, w, a)
    mu2hat <- sdev^2 * mu2hat
    jvhat.d11  <- sum(mu2hat)
     
    mutilde <- postmed(x, w, prior = prior, a = a)
    mutilde <- sdev * mutilde
    jvtilde.js05 <- sum(mutilde^2)
  }
 
  # my additions. pj.m, pj.b (posterior jump probabilities)
  if (pr == "l") {
    sx <- sign(x)
	  x <- abs(x)
    wpost <- 1 - (1 - w)/(1 + w * beta.laplace(x, a))
    ef <- exp(pmin(2 * a * x, 100))
    pj.m <- wpost / (1 + ef * pnorm(-x - a) / pnorm(x - a)) 
    pj.b <- wpost
    x <- sx * x  # get back the original x. 
    jvtilde.d11 <- jv.median(x, a, wpost, sdev)
  }
 
	if(!verbose)
		return(muhat)
	retlist <- list(muhat = muhat, x = x, threshold.sdevscale = tt, 
		threshold.origscale = tcor, prior = prior, w = w, a = a, wpost = wpost,
		bayesfac = bayesfac, sdev = sdev, threshrule = threshrule,
    mu2hat = mu2hat, jvhat.d11 = jvhat.d11, muhat2 = muhat2, 
    jvhat.js05 = jvhat.js05, pj.b = pj.b, pj.m = pj.m, jvtilde.d11 = jvtilde.d11, 
    jvtilde.js05 = jvtilde.js05, mutilde = mutilde)
	if(pr == "c")
		retlist <- retlist[-7]
	if(threshrule == "none")
		retlist <- retlist[-1]
	return(retlist)
}
