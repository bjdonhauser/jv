"ret.emp.lm08" <- function(lret.mat) {  
  #
  #
  
  # Initializing local parameters
  n.obs       <- dim(lret.mat)[1]  # 389 for 1min returns
  n.days      <- dim(lret.mat)[2]  # 60 for 2008.01--2008.03 
  
  model.names <- c("lm08.l.05", "lm08.l.01", 
                   "lm08.g.05", "lm08.g.01", 
                   "lm08.g.naive")
  n.models    <- length(model.names)
  
  # Initialzing return variables
  jvhat       <- matrix(rep(0, n.days * n.models), 
                        nrow = n.days, ncol = n.models,
                        dimnames = list(NULL, model.names))                    
  ivhat       <- matrix(rep(0, n.days * n.models), 
                        nrow = n.days, ncol = n.models,
                        dimnames = list(NULL, model.names))
  
  # Calculate Realized Measures
  rv <- colSums(lret.mat^2)
  
  for (i in 1:n.days) {
    
    ret.list <- lm08(lret.mat[ , i], stat.model = "local", 
                     alpha = 0.05, verbose = T)
    jvhat[i, model.names[1]] <- ret.list$jvhat
    
    ret.list <- lm08(lret.mat[ , i], stat.model = "local", 
                     alpha = 0.01, verbose = T)
    jvhat[i, model.names[2]] <- ret.list$jvhat
    
    ret.list <- lm08(lret.mat[ , i], stat.model = "global", 
                     alpha = 0.05, verbose = T)
    jvhat[i, model.names[3]] <- ret.list$jvhat
    
    ret.list <- lm08(lret.mat[ , i], stat.model = "global", 
                     alpha = 0.01, verbose = T)
    jvhat[i, model.names[4]] <- ret.list$jvhat
    
    jvhat[i, model.names[5]] <- lm08.naive(lret.mat[ , i], stat.model="global")
    
    ivhat[i, ] <- rv[i] - jvhat[i, ]
  }
  
  model.ret <- list(rv = rv, jvhat = jvhat, ivhat = ivhat)
  return(model.ret)
}
