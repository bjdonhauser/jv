"lm08.b2a" <- function(beta, stat.model = "global"){
  alpha <- NA
  if (stat.model == "global"){
    alpha <- -exp(-exp(-beta)) + 1   
  } else if (stat.model == "local"){
    alpha <- (1 - pnorm(abs(beta))) * 2    
  }
  return(alpha)
}