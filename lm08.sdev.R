"lm08.sdev" <- function(x, sdev = NA, 
                        window.type = "full", window.len = NA) {
  # Scale x to have unit variance 
  #
  # Args: 
  #   x: Data. To be interpreted as log-returns.
  #      Implies preprocessing on log-prices. 
  #   sdev: [NA | "sv", "rbpv", "mad", numeric]
  #   window.type: ["smoothed" | "predictive", "filtered"]
  #   window.len: [NA]
  #
  # Returns:
  #   sdev: estimate of the standard deviation of an observation
  #
  if (window.type != "full" && sdev != "sv") {
    error("Error: cannot have `window.type` non-`full` without `sdev=sv`")
  }
  if (sdev == "mad" || is.na(sdev)) {
    sdev <- mad(x, center = 0)
  } else if (sdev == "rbpv") {
    n <- length(x)
    c <- sqrt(2/pi)
    sdev <- (1/c) * sqrt(mean(abs(x[1:(n-1)]*x[2:n])))  # note the `mean()`!
  } else if (sdev == "sv") {
    sdev <- lm08.sdev.sv(x, window.type, window.len)
  }
  # Otherwise it's assumed that `sdev` passed is numeric and returns itself
  # This can be interpreted as having a prior on `sdev`
  return(sdev)
}