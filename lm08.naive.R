"lm08.naive" <- function(x, stat.model = "global", sdev = "rbpv") {
  # Returns my naive shrinkage estimator of jump variation
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
  #
  # Returns:
  #   jv.naive
  #
  # Requires:
  #   EbayesThresh, lm08.tfromx(), lm08.sdev()
  sdev <- lm08.sdev(x, sdev)
  x.unit <- x / sdev
  n <- length(x)
  c <- sqrt(2/pi)
  s.n <- 1/(c*sqrt(2*log(n)))  
    # same as lm08
  c.n <- (s.n * 2 * log(n)) - ((s.n/2) * (log(pi) + log(log(n))))  
    # same as lm08
  x.stat <- (abs(x.unit) - c * c.n) / (c * s.n)  
    # `*c` because of scaling in x.unit
  x.weight <- exp(-exp(-x.stat))
  jv.naive <- sum((x.weight * x)^2)
  return(jv.naive)
}