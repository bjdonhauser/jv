"lm08.thresh" <- function(alpha, stat.model, n = NA){
  # See lm08.tfromx for basically the same function
  if (stat.model == "local") {
    t <- qnorm(1-alpha/2)  
    # 2.58 for alpha = 0.01. 
    # Value n irrelevant because this is a local procedure 
  } else if (stat.model == "global") {
    if (is.na(n)) {
      stop("Error: If stat.model == global, then must pass a value for n")
    }
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