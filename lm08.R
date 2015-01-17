"lm08" <- function(x, stat.model = "global", alpha = 0.01, sdev = "rbpv", 
                   threshrule = "hard", verbose = FALSE) {
  # Performs the analysis of LM08
  #
  # Args: 
  #   x: Data. To be interpreted as log-returns.
  #      Implies preprocessing on log-prices. 
  #      Also, have been scaled by the estimate of sdev!
  #   stat.model: ["global"|"local"] specify whether the threshold returned
  #               is local or includes the global, max order statistic 
  #               correction. 
  #   alpha: [0.01] the significance level to which the threshold corresponds
  #   sdev: ["rbpv" | "mad", numeric]
  #   threshrule: ["hard" | "soft", "none"]
  #   verbose: [FALSE | TRUE]
  #
  # Returns:
  #   muhat: estimation of jumps
  #   ret.list: the full list of output, including muhat
  #
  # Requires:
  #   EbayesThresh, lm08.tfromx(), lm08.sdev()
  #
  sdev <- lm08.sdev(x, sdev)
  x.unit <- x / sdev
  t.unit <- lm08.tfromx(x.unit, stat.model, alpha)
  jumps.loc <- abs(x.unit) > t.unit
  if (threshrule == "hard") {
    muhat.unit <- threshld(x.unit, t.unit, hard = TRUE)
  } else if (threshrule == "soft") {
    muhat.unit <- threshld(x.unit, t.unit, hard = FALSE)
  } else if (threshrule == "none") {
    muhat.unit <- NA
  }
  muhat <- muhat.unit * sdev
  if(!verbose) {
    return(muhat)
  } else if (verbose) {
    w <- wfromt(t.unit)
    t <- t.unit * sdev
    jvhat <- sum(muhat[jumps.loc]^2)
    ret.list <- list(muhat = muhat, muhat.unit = muhat.unit, 
                     x = x, x.unit = x.unit, 
                     t = t, t.unit = t.unit, 
                     sdev = sdev, w = w, jvhat = jvhat,
                     jumps.loc = jumps.loc, threshrule = threshrule)
    return(ret.list)
  }
}