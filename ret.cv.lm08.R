"ret.cv.lm08" <- function(n.sim = 1000,   # 1000, 5000
                          n.jumps = 0,    # 0, 3, 10, 30                    
                          jump.scale = 7  # 4, 7, 15
                          ) {  
  #
  set.seed(1)
  n.obs      <- 390
  vol        <- .20         # .20, .30, .60
  delta      <- 1/(390*252)
  
  iv <- rep(0, n.sim)
  jv <- rep(0, n.sim)
  qv <- rep(0, n.sim)
  rv <- rep(0, n.sim)
  
  model.names <- c("lm08.l.05", "lm08.l.01", 
                   "lm08.g.05", "lm08.g.01", 
                   "lm08.g.naive", "bns04") 
  n.models    <- length(model.names)  
  jvhat       <- matrix(rep(0, n.sim * n.models), 
                        nrow = n.sim, ncol = n.models,
                        dimnames = list(NULL, model.names))                    
  ivhat       <- matrix(rep(0, n.sim * n.models), 
                        nrow = n.sim, ncol = n.models,
                        dimnames = list(NULL, model.names))
  jvhat.err   <- matrix(rep(0, n.sim * n.models), 
                        nrow = n.sim, ncol = n.models,
                        dimnames = list(NULL, model.names))                      
  ivhat.err   <- matrix(rep(0, n.sim * n.models), 
                        nrow = n.sim, ncol = n.models,
                        dimnames = list(NULL, model.names))
  
  for (i in 1:n.sim) {
    eps <- rnorm(n = n.obs, sd = vol * sqrt(delta))
    mu  <- jump.scale * vol * sqrt(delta) * 
      sample(c(runif(n.jumps, -1, 1), rep(0, n.obs - n.jumps)))
    x   <- mu + eps
    
    iv[i] <- vol^2 * delta * n.obs
    jv[i] <- sum(mu^2)
    qv[i] <- iv[i] + jv[i]
    rv[i] <- sum(x^2)  
    
    # jv
    ret.list <- lm08(x, stat.model = "local", alpha = 0.05, verbose = T)
    jvhat[i, model.names[1]] <- ret.list$jvhat
    ret.list <- lm08(x, stat.model = "local", alpha = 0.01, verbose = T)
    jvhat[i, model.names[2]] <- ret.list$jvhat
    ret.list <- lm08(x, stat.model = "global", alpha = 0.05, verbose = T)
    jvhat[i, model.names[3]] <- ret.list$jvhat
    ret.list <- lm08(x, stat.model = "global", alpha = 0.01, verbose = T)
    jvhat[i, model.names[4]] <- ret.list$jvhat
    jvhat[i, model.names[5]] <- lm08.naive(x, stat.model="global")

    
    # non-bns04: jvhat.err, ivhat, ivhat.err                 
    jvhat.err[i, model.names[-n.models]] <- jvhat[i, model.names[-n.models]] - jv[i] 
    ivhat[i, model.names[-n.models]]     <- rv[i] - jvhat[i, model.names[-n.models]]
    ivhat.err[i, model.names[-n.models]] <- ivhat[i, model.names[-n.models]] - iv[i]
    
    # bns04: jvhat, jvhat.err, ivhat. ivhat.err
    ret.list <- bns04(x)
    jvhat[i, model.names[n.models]] <- ret.list$jvhat    
    jvhat.err[i, model.names[n.models]] <- jvhat[i, model.names[n.models]] - 
                                             jv[i]
    ivhat[i, model.names[n.models]]     <- ret.list$ivhat
    ivhat.err[i, model.names[n.models]] <- ivhat[i, model.names[n.models]] - 
                                             iv[i]
    
  }
  
  model.ret <- list(iv = iv, jv = jv, qv = qv, rv = rv,
                    jvhat = jvhat, ivhat = ivhat, 
                    jvhat.err = jvhat.err, ivhat.err = ivhat.err)
  
  #################### Cleanup
  rm(i, n.obs, vol, delta, jump.scale, n.jumps, n.sim, 
     iv, jv, qv, rv,
     model.names,
     jvhat, ivhat, jvhat.err, ivhat.err,
     eps, mu, x,
     ret.list)
  
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