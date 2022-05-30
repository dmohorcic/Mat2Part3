HMC <- function(f, df, M, step.size, L, q.0, m) {
  # f is proportional to our distribution
  # df is the gradient of f
  # M is (diagonal) mass matrix
  # step.size is the step size
  # L is the number of steps
  # q.0 is the target density
  # m is the number of samples we want
  
  n <- length(q.0)
  q.samples <- array(dim = c(m, n)) # store all samples from f
  
  # initialize the system
  M.inv = solve(M)
  q.prev <- q.0
  
  for(i in 1:m) {
    # sample new momentum
    p <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = M)
    
    # init previous value and current position
    p.prev <- p
    q <- q.prev
    
    # leapfrog
    p <- p - step.size/2 * df(q)
    for(l in 1:L) {
      q <- q + step.size * p
      if(l != L) {p <- p - step.size * df(q)}
    }
    p <- p - step.size/2 * df(q)
    p <- -p
    
    # calculate acceptance probabilities
    H.start <- as.numeric(0.5*(p.prev%*%M.inv%*%p.prev) + f(q.prev))
    H.end <- as.numeric(0.5*(p%*%M.inv%*%p) + f(q))
    alpha <- min(1, exp(H.start - H.end))
    u <- runif(1, 0, 1)
    
    # save sample
    if(u <= alpha) {
      q.prev <- q
    }
    q.samples[i,] <- q.prev
  }
  
  return(q.samples)
}
