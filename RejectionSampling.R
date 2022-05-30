RS <- function(f, M, a, b, n, m) {
  # f is proportional to our distribution
  # M is constant so that M >= max(f)
  # a, b represent the range from which to sample x~Uniform(a, b)
  #      both can be numbers or vectors (of size n)
  # n is the number of dimensions
  # m is the number of samples we want
  
  if(length(a) == 1) {a <- rep(a, n)}
  if(length(b) == 1) {b <- rep(b, n)}
  
  x.samples <- array(dim=c(m, n))
  
  # keep track of number of found samples
  m.found <- 0
  
  # iterate until enough samples are found
  while(m.found < m) {
    x.prop <- runif(n, a, b) # use uniform wrapper for proposed sample
    u <- runif(1, 0, 1)
    if(u <= prod(b-a)/M * f(x.prop)) { # accept the sample
      m.found <- m.found + 1
      x.samples[m.found,] <- x.prop
    }
  }
  
  return(x.samples)
}