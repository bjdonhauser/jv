"ret.sv" <- function(jump.scale = 7,
                     n.sim = 1000, 
                     n.jumps = 0, 
                     sdev.parm = "mad", 
                     window.type = "full", 
                     window.len = NA) {
  #
  set.seed(1)
  n.obs      <- 390
  vol        <- .20  # 19.774%, 20%, 30%, 60% ann. 
                     # => LM08 vals, 6.38e-4, 9.57e-4, 1.91e-3 
  delta      <- 1/(390*252) 

  theta    <- rep(0, 3)
  theta[2] <- 8.0369                        # Reversion parm. Set to LM08 sim val
  theta[1] <- theta[2] * vol^2              # Set to target LRVar
  theta[3] <- sqrt((2/ 3.4390) * theta[1])  # Set so that df of marginal constant

  #################### Simulate
  ret <- lm08.sim.sv(N = n.obs, M = n.sim, delta = delta, theta = theta)

  mu  <- matrix(rep(0, n.obs * n.sim), nrow = n.obs, ncol = n.sim)
  for (i in 1:n.sim) {
    mu[ , i]  <- jump.scale * vol * sqrt(delta) * 
                 sample(c(runif(n.jumps, -1, 1), rep(0, n.obs - n.jumps)))
  }      
  
  x   <- mu + ret$logr

  #################### Calculate Realized Measures
  iv <- delta * colSums(ret$sig2)
  jv <- colSums(mu^2)
  qv <- iv + jv
  rv <- colSums(x^2)
  
  #################### Do EBT Analysis
  model.names <- c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11", "bns04")  

  jvhat       <- matrix(rep(0, n.sim * length(model.names)), 
                        nrow = n.sim, ncol = length(model.names),
                        dimnames = list(NULL, model.names))                    
  ivhat       <- matrix(rep(0, n.sim * length(model.names)), 
                        nrow = n.sim, ncol = length(model.names),
                        dimnames = list(NULL, model.names))
  jvhat.err   <- matrix(rep(0, n.sim * length(model.names)), 
                        nrow = n.sim, ncol = length(model.names),
                        dimnames = list(NULL, model.names))                      
  ivhat.err   <- matrix(rep(0, n.sim * length(model.names)), 
                        nrow = n.sim, ncol = length(model.names),
                        dimnames = list(NULL, model.names))

  for (i in 1:n.sim) {
    ret.list <- ebayesthresh.mod(x[ , i], threshrule = "mean", verbose = T, 
                                 sdev = sdev.parm, 
                                 window.type = window.type, 
                                 window.len = window.len)
    jvhat[i, c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11")] <- 
        unlist(ret.list[c("jvhat.js05", "jvtilde.js05", "jvhat.d11", 
                        "jvtilde.d11")])
    jvhat.err[i, c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11")] <- 
      jvhat[i, c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11")] - jv[i] 
    ivhat[i, c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11")] <- 
      rv[i] - jvhat[i, c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11")]
    ivhat.err[i, c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11")] <- 
      ivhat[i, c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11")] - iv[i]
                       
    ret.list <- bns04(x[ , i])
    jvhat[i, "bns04"]     <- ret.list$jvhat
    jvhat.err[i, "bns04"] <- jvhat[i, "bns04"] - jv[i]
    ivhat[i, "bns04"]     <- ret.list$ivhat
    ivhat.err[i, "bns04"] <- ivhat[i, "bns04"] - iv[i]
  
    #cat(i, ", ", sep = "")
  }

  model.ret <- list(iv = iv, jv = jv, qv = qv, rv = rv,
                    jvhat = jvhat, ivhat = ivhat, 
                    jvhat.err = jvhat.err, ivhat.err = ivhat.err)

  #################### Cleanup
  rm(n.obs, vol, delta, jump.scale, n.jumps, n.sim, 
     theta, ret, mu, x,
     iv, jv, qv, rv,
     model.names, ret.list,
     jvhat, ivhat, jvhat.err, ivhat.err)

  #################### Interp
#===============================#
# jvhat error                 
#===============================#
# ME
model.ret$jvhat.me <- matrix(colMeans(model.ret$jvhat.err), ncol=1, 
       dimnames = list(colnames(model.ret$jvhat.err), NULL))
# MAE
model.ret$jvhat.mae <-matrix(colMeans(abs(model.ret$jvhat.err)), ncol=1, 
       dimnames = list(colnames(model.ret$jvhat.err), NULL))
# RMSE
model.ret$jvhat.rmse <-matrix(colMeans(sqrt(model.ret$jvhat.err^2)), ncol=1, 
       dimnames = list(colnames(model.ret$jvhat.err), NULL))
       
#===============================#
# jvhat error as % of qv                 
#===============================#
# MPE
model.ret$jvhat.mpe  <-matrix(colMeans(model.ret$jvhat.err/model.ret$qv), ncol=1, 
       dimnames = list(colnames(model.ret$jvhat.err), NULL))
# MAPE
model.ret$jvhat.mape <-matrix(colMeans(abs(model.ret$jvhat.err/model.ret$qv)), ncol=1, 
       dimnames = list(colnames(model.ret$jvhat.err), NULL))
# RMSPE
model.ret$jvhat.rmspe <-matrix(sqrt(colMeans((model.ret$jvhat.err/model.ret$qv)^2)), ncol=1, 
       dimnames = list(colnames(model.ret$jvhat.err), NULL))

#===============================#
# ivhat error                 
#===============================#
# ME
model.ret$ivhat.me <- matrix(colMeans(model.ret$ivhat.err), ncol=1, 
       dimnames = list(colnames(model.ret$ivhat.err), NULL))
# MAE
model.ret$ivhat.mae <-matrix(colMeans(abs(model.ret$ivhat.err)), ncol=1, 
       dimnames = list(colnames(model.ret$ivhat.err), NULL))
# RMSE
model.ret$ivhat.rmse <-matrix(sqrt(colMeans(model.ret$ivhat.err^2)), ncol=1, 
       dimnames = list(colnames(model.ret$ivhat.err), NULL))
       
#===============================#
# ivhat error as % of qv                 
#===============================#
# MPE
model.ret$ivhat.mpe  <- matrix(colMeans(model.ret$ivhat.err/model.ret$qv), 
                               ncol=1, 
                               dimnames = list(colnames(model.ret$ivhat.err), 
                                               NULL))
# MAPE
model.ret$ivhat.mape <- matrix(colMeans(abs(model.ret$ivhat.err/model.ret$qv)), 
                               ncol=1, 
                               dimnames = list(colnames(model.ret$ivhat.err), 
                                               NULL))
# RMSPE
model.ret$ivhat.rmspe <- matrix(sqrt(colMeans((model.ret$ivhat.err/model.ret$qv)^2)), 
                                ncol=1, 
                                dimnames = list(colnames(model.ret$ivhat.err), 
                                                NULL))
return(model.ret)
}