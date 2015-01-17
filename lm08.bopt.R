"lm08.bopt" <- function(n = 390, N = 3, sigma = 1, T = 1){
  # See LeeMyk08 Thm5
  beta.opt <- -log(sigma * sqrt(T) * N / sqrt(2 * n * log(n)))
  return(beta.opt)
}