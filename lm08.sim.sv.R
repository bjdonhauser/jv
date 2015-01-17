"lm08.sim.sv" <- function(N = 390, M = 1, n.burn = 200, delta = 1/98280, 
                          theta = c(0.3141, 8.0369, 0.43)) {
  # Simulate log-returns and variance from the CIR SV-model used in LM08
  # The LRVol implied by the parameterization is 0.013day, 0.20year
  #
  # Args: 
  #   N: [390] number of simulation steps. Default of 390 = 6.5hrs*60min =
  #      the number of minutes in a trading day
  #   M: [1] number of trajectories. For a typical test would use M=1000
  #   n.burn: [200] number of burn-in steps. 
  #   delta: [1/98280] Proportion of a full step taken. Since the parameters
  #          are annualized and there are 252days*6.5hrs*60min = 98280 minutes
  #          in a trading year. The corresponding data respond to minute steps
  #   theta: [c(0.3141, 8.0369, 0.43))] 
  #          Defaults correspond to c(theta_1, theta_2, theta_3) in Iac08. 
  #          LRVol = sqrt(theta_1/theta_2) 
  #                = sqrt(0.04) = 0.20 (annualized) in the default case. 
  #          theta_1 = 0.3141, 0.7233, 2.8933 correspond to LRVol
  #                    0.20,   0.30,     0.60,    when theta_2 = 8.0369
  #          theta_3 = 0.7925, 1.2027, 2.4055, corresponding to the theta1 vals
  #                    are max values for stat. cond.: 2*theta1 > theta3^2
  #                    Default: 2*theta1/theta3^2 = 2*0.3141/0.43^2 = 3.4 
  #                    Holding this ratio constant effectively keeps the df
  #                    in the non-central chi-squared marginal dbtn the same
  #                    and gives:
  #          theta_3 = 0.43,  0.6523, 1.3046
  #          
  #
  # Returns:
  #   sim
  #     sim$sig2: NxM matrix of variances in annualized terms
  #     sim$logr: NxM matrix of log-returns
  # 
  # Uses:
  #   sde
  #
  #
  sig2 <- sde.sim(N = (N + n.burn), M = M, delta = delta, 
                  X0 = theta[1]/theta[2], theta = theta, model = "CIR")
  # theta[1]/theta[2] corresponds to the long-run value of the variance
  sig2 <- as.matrix(sig2)[(n.burn + 2):(N + n.burn + 1), 1:M]  # NxM var matrix
  logr <- rnorm(n = N * M, sd = sqrt(sig2 * delta))
  logr <- matrix(logr, nrow = N, ncol = M)
  #logp <- cumsum(c(log(100), logr))  # Run if log-price series is needed
  ret.list <- list(logr = logr, sig2 = sig2)
  return(ret.list)
}