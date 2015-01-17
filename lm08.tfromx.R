"lm08.tfromx" <- function(x, stat.model = "global", alpha = 0.01) {
  # Calculate the threshold above which (in absolute value) x is taken to have
  # a jump embedded. 
  # This function is itended to be called on x with variance 1. 
  # So it has been pre-scaled by an estimate of the LR stdev of x. 
  # That can be a mad-based estimate or the rbpv-based estimate as in LM08.
  # Note that these thresholds don't really require the x-values. 
  # That's because the x passed has already been stabilized to unit variance. 
  #
  # Args: 
  #   x: Data. To be interpreted as log-returns.
  #      Implies preprocessing on log-prices. 
  #      Also, have been scaled by the estimate of sdev!
  #   stat.model: ["global"|"local"] specify whether the threshold returned
  #               is local or includes the global, max order statistic 
  #               correction. 
  #   alpha: [0.01] the significance level to which the threshold corresponds
  #
  # Returns:
  #   t: a positive scalar giving the threshold above which scaled x in
  #      absolute value is consider to have a jump. 
  #
  n <- length(x)
  if (stat.model == "local") {
    t <- qnorm(1-alpha/2)  
    # 2.58 for alpha = 0.01. 
    # Value n irrelevant because this is a local procedure 
  } else if (stat.model == "global") {
      c <- sqrt(2/pi)
      s.n <- 1/(c*sqrt(2*log(n)))  
      c.n <- (s.n * 2 * log(n)) - ((s.n/2) * (log(pi) + log(log(n))))
      beta <- -log(-log(1-alpha))  # 4.60 for alpha = 0.01
      t <- c * (c.n + s.n * beta)  # 4.36 for n = 390, alpha = 0.01
                                   # 3.48 if we'd forgotten to scale b `c`
  } else {
      stop("Error: passed bad value for `stat.model`")
  }
  return(t)  
}