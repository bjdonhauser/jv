"ret.emp.bns04" <- function(lret.mat) {
  #
  #
  
  # Initializing local parameters
  n.obs       <- dim(lret.mat)[1]  # 389 for 1min returns
  n.days      <- dim(lret.mat)[2]  # 60 for 2008.01--2008.03
  
  # Initialzing return variables
  jvhat       <- matrix(rep(0, n.days), nrow = n.days, 
                        dimnames = list(NULL, "bns04"))                    
  ivhat       <- matrix(rep(0, n.days), nrow = n.days,
                        dimnames = list(NULL, "bns04"))
  # Calculate Realized Measures
  rv <- colSums(lret.mat^2)
  
  for (i in 1:n.days) {
    ret.list <- bns04(wmt.lret.mat[ , i])
    jvhat[i, "bns04"] <- ret.list$jvhat
    ivhat[i, "bns04"] <- ret.list$ivhat
  }
  
  model.ret <- list(rv = rv, jvhat = jvhat, ivhat = ivhat)
  return(model.ret)
}