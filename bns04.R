"bns04" <- function(x) {
  # Performs analysis of bns04
  #
  # Args: 
  #   x: Data. To be interpreted as log-returns.
  #      Implies preprocessing on log-prices.
  #
  # Returns:
  #   ret: 
  #   ret$rv
  #   ret$ivhat
  #   ret$jvhat
  #
  c <- sqrt(2/pi)
  n <- length(x)
  rv <- sum(x^2)
  ivhat <- (1/c^2)*sum(abs(x[2:n])*abs(x[1:(n-1)]))
  jvhat <- rv - ivhat
  ret.list <- list(rv = rv, ivhat = ivhat, jvhat = jvhat)
}
