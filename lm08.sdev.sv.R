"lm08.sdev.sv" <- function(x, window.type = "full", window.len = NA) {
  # Estimate sdev by a local, windowed version of the MAD estimator
  #
  # Args: 
  #   x: Data. To be interpreted as log-returns.
  #      Implies preprocessing on log-prices. 
  #   window.type: ["smoothed" | "filtered", "full"]
  #   window.len: [NA] Normal values: 5, 10, 15, 30 (for minutes)
  #
  # Returns:
  #   sdev: estimate of the standard deviation of an observation
  #
  if ((window.type != "full") && is.na(window.len)) 
    stop("Error: Must supply a window length `window.len` if not using the 
         `full` window")
  n.obs <- length(x)
  sdev  <- rep(0, n.obs) 
  if (window.type == "full") {
    # Uses entire day's data to estimate one spot volatility for the entire day
    sdev <- rep(mad(x, center = 0), n.obs)
  } else if (window.type == "filtered") {
    for (i in 1:n.obs) {
      j <- max(1, i - (window.len - 1))
      sdev[i] <- mad(x[j:i], center=0)
    }
  } else if (window.type == "smoothed") {
    for (i in 1:n.obs) {
      j <- max(1, i - (window.len - 1))
      k <- min(n.obs, i + (window.len - 1))
      sdev[i] <- mad(x[j:k], center=0)
    }      
  } else {
      stop("Error: passed a bad value for `window.type` parameter.
            User is trying to indicate a window.type that does not exist.")
  }
  return(sdev)
}