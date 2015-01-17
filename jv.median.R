"jv.median" <- function(x, a, wpost, sdev, n.sample = 200) {
  # Get the median estimate of JV
  #
  # Args: 
  #   x: Data. To be interpreted as log-returns. 
  #      Implies preprocessing on log-prices. 
  #   a: estimated or prior a in the fitted Laplace distribution
  #   wpost: posterior estimate of jump probability (same as pj.b)
  #   sdev: estimate of standard deviation. may be a scalar or n.obsx1 vector
  #         in the case of stochastic volatility estimation.
  #   n.sample: [200] number of samples from the individual posterior 
  #             distributions. 
  #             200 <=> 0.1% changes in median estimate
  #             3000 <=> 0.01% changes in median estimate
  #             
  # Returns:
  #   ret: 
  #
  # Requires:
  #   packages: rtmvnorm
  #
  n.obs <- length(x)
  sample <- matrix(rep(0, n.obs * n.sample), nrow = n.obs, ncol = n.sample)
  cp1 <- exp(-a * x) * pnorm(x - a)
  cp2 <- exp(a * x) * pnorm(-(x + a))
  c <- cp1 + cp2
  pj0 <- 1 - wpost  # probability of no jump
  pj1 <- wpost * cp1 / c
  pj2 <- wpost * cp2 / c
  for (i in 1:n.obs) {  
    n.bucket <- rmultinom(1, n.sample, c(pj0[i], pj1[i], pj2[i]))
    if (n.bucket[2] > 0) {
      tn1 <- rtmvnorm(n.bucket[2], mean = (x[i] - a), sigma = 1, lower = 0,
                      algorithm = "gibbs")
    } else {
      tn1 <- integer(0)
    }
    if (n.bucket[3] > 0) {
      tn2 <- rtmvnorm(n.bucket[3], mean = (x[i] + a), sigma = 1, upper = 0, 
                      algorithm = "gibbs")  
    } else {
      tn2 <- integer(0)
    }
    sample[i, ] <- sample(c(rep(0, n.bucket[1]), tn1, tn2))
  }
  jv.sample <- colSums((sdev * sample)^2)  # (n.sample x 1) vector
  jv.median <- median(jv.sample)
}

