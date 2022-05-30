MH <- function(f, Q, m, x.0=NULL) {
  # f is proportional to our distribution
  # Q is the transition density function
  # m is the number of samples we want
  # x.0 is the starting point (unused)
  
  n <- sqrt(length(Q))
  x.samples <- array(dim=c(m, n)) # store all samples from f
  
  # Multivariate normal proposal for starting point
  if(!is.null(x.0)) {
    x.prev <- x.0
  } else {
    x.prev <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = Q)
  }
  
  for(i in 1:m) {
    # sample candidate state
    x.prop <- MASS::mvrnorm(1, mu = x.prev, Sigma = Q)
    
    # accept candidate with some probability
    alpha <- min(1, f(x.prop)/f(x.prev))
    
    u <- runif(1, 0, 1)
    if(u <= alpha) {
      x.prev <- x.prop
    }
    
    # save sample
    x.samples[i,] <- x.prev
  }
  
  return(x.samples)
}
