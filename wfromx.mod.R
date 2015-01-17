"wfromx.mod" <-
function(x, prior = "laplace", a = 0.5)
{
#  given the vector of data x and the function betaf
#   which finds beta(x), 
#  find the value of w that zeroes S(w) in the
#  range 
#
#  works by successive bisection, carrying out nits harmonic bisections
#   of the original interval between wlo and 1.  
#  
	pr <- substring(prior, 1, 1)
	tuniv <- sqrt(2 * log(length(x)))
  #wlo <- wfromt(tuniv, prior, a)          # standard
  wlo <- 0                                # mod
	#wlo <- wfromt(tuniv, prior, a) * 1e-14  # the modification!
	if(pr == "l") {
		beta <- beta.laplace(x, a)
	}
	if(pr == "c") {
		beta <- beta.cauchy(x)
	}
	whi <- 1
	beta <- pmin(beta, 1e20) 
	shi <- sum(beta/(1 + beta))
	if(shi >= 0)
		return(w = 1)
	slo <- sum(beta/(1 + wlo * beta))
	if(slo <= 0)
		return(w = wlo)
 
  wlo <- wfromt(tuniv, prior, a) * 1e-14  # MOD. Consider changing to 1e-8
 
	for(j in (1:30)) {
    #cat("j: ", j, "\n")
    #cat("wlo: ", wlo, "\n")
    #cat("whi: ", whi, "\n")
		wmid <- sqrt(wlo * whi)
		smid <- sum(beta/(1 + wmid * beta))
		if(smid == 0)
			return(w = wmid)
		if(smid > 0) {
			wlo <- wmid
		}
		else {
			whi <- wmid
		}
	}
	return(w = sqrt(wlo * whi))
}
