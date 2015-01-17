"lm08.a2b" <- function(alpha, stat.model = "global"){
  beta <- NA
  if (stat.model == "global"){
    beta <- -log(-log(1-alpha))    
  } else if (stat.model == "local"){
    beta <- qnorm(1-alpha/2)    
  }
  return(beta)
}