"ret.emp.ebt" <- function(lret.mat, sdev.parm = "mad", 
                          window.type = "full", window.len = NA) {  
  #
  #
  
  # Initializing local parameters
  n.obs       <- dim(lret.mat)[1]  # 389 for 1min returns
  n.days      <- dim(lret.mat)[2]  # 60 for 2008.01--2008.03 
  
  model.names <- c("hat.js05", "tilde.js05", "hat.d11", "tilde.d11")
  n.models    <- length(model.names)
  
  # Initialzing return variables
  jvhat       <- matrix(rep(0, n.days * n.models), 
                        nrow = n.days, ncol = n.models,
                        dimnames = list(NULL, model.names))                    
  ivhat       <- matrix(rep(0, n.days * n.models), 
                        nrow = n.days, ncol = n.models,
                        dimnames = list(NULL, model.names))
  sdev        <- matrix(rep(0, n.obs * n.days), nrow = n.obs, ncol = n.days)
  muhat       <- matrix(rep(0, n.obs * n.days), nrow = n.obs, ncol = n.days)
  mutilde     <- matrix(rep(0, n.obs * n.days), nrow = n.obs, ncol = n.days)
  
  # Calculate Realized Measures
  rv <- colSums(lret.mat^2)
  
  for (i in 1:n.days) {
    ret.list <- ebayesthresh.mod(wmt.lret.mat[ , i], 
                                 threshrule = "mean", 
                                 sdev = sdev.parm, 
                                 window.type = window.type, 
                                 window.len = window.len,
                                 verbose = T)
    
    jvhat[i, model.names] <- unlist(ret.list[c("jvhat.js05", "jvtilde.js05", 
                                               "jvhat.d11", "jvtilde.d11")])
    ivhat[i, model.names] <- rv[i] - jvhat[i, model.names]
    
    sdev[ , i]    <- ret.list$sdev
    muhat[ , i]   <- ret.list$muhat
    mutilde[ , i] <- ret.list$mutilde
    
    cat(i, ", ", sep = "")
  }
  
  model.ret <- list(rv = rv, jvhat = jvhat, ivhat = ivhat, 
                    sdev = sdev, muhat = muhat, mutilde = mutilde)
  return(model.ret)
}
