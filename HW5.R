library(mcmcse)
library(ggplot2)


###################
# Implementations #
###################

source("MetropolisHastings.R")
source("HamiltonianMC.R")
source("RejectionSampling.R")

plot.samples <- function(samples) {
  ggplot(data=samples, aes(x=x1, y=x2, color=chain))+
    geom_path(alpha=0.3)+geom_point(size=1)
}

plot.samples.x1 <- function(samples) {
  c <- length(unique(samples$chain))
  t <- rep(1:(nrow(samples)/c), c)
  ggplot(data=samples, aes(x=t, y=x1, color=chain))+
    geom_line(alpha=0.3)+geom_point(size=1)
  
}

plot.samples.x2 <- function(samples) {
  c <- length(unique(samples$chain))
  t <- rep(1:(nrow(samples)/c), c)
  ggplot(data=samples, aes(x=t, y=x2, color=chain))+
    geom_line(alpha=0.3)+geom_point(size=1)
}


########
# Code #
########

chains <- 5
m <- 1000


################################
# a) bivariate standard normal #
################################

Q = matrix(data=c(1, 0, 0, 1), ncol=2)
Q.INV <- solve(Q)
x.0 = c(0, 0)

f <- function(x) {
  res <- 1/(2*pi) * 1/sqrt(det(Q)) * exp(-0.5*(x-x.0)%*%Q.INV%*%(x-x.0))
  return(as.numeric(res))
}
df <- function(x) {
  val <- f(x)
  res <- -val * t(Q.INV%*%x)
  return(as.numeric(res))
}


# with Metropolis-Hastings

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- MH(f, Q, m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}


plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)


# with Hamiltonian MC

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- HMC(f, df, M = Q, step.size = 0.012, L = 25, q.0 = c(0, 0), m = m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}

plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)


ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)


# grid for best params
best.ess <- 0
best.step <- 0
best.L <- 0

for(step.size in seq(0.002, 0.01, 0.002)) {
  for(L in seq(10, 50, 10)) {
    cat(step.size, L, "\n")
    # sample
    samples <- NULL
    times <- rep(0, chains)
    for(c in 1:chains) {
      start <- as.numeric(Sys.time())
      sample <- HMC(f, df, M = Q, step.size = step.size, L = L, q.0 = c(0, 0), m = m)
      end <- as.numeric(Sys.time())
      
      samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
      times[c] <- end-start
    }
    
    ESS <- array(dim = c(chains, 2))
    ESS.s <- array(dim = c(chains, 2))
    for(c in 1:chains) {
      chain.n <- samples[samples$chain == as.character(c),][,1:2]
      ESS[c,] <- ess(chain.n) # ESS
      ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
    }
    
    if(mean(ESS) > best.ess) {
      best.ess <- mean(ESS)
      best.step <- step.size
      best.L <- L
    }
  }
}


# with rejection sampling

a <- rep(-5, 2)
b <- rep(5, 2)
f.max <- dnorm(0)^2 # we know it for bivariate standard normal
M <- f.max*prod(b-a)

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- RS(f, M, a, b, 2, m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}

plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)


######################
# b) banana function #
######################

B <- 0.05
minus_logf <- function(x) {
  -(-(x[1]^2)/200 - 0.5*(x[2] + B*x[1]^2 - 100*B)^2)
}
minus_logf_grad <- function(x) {
  g1 <- -(x[1])/100- 1.0 * (2* B * x[1]) * (x[2]+ B * x[1]^2 - 100*B)
  g2 <- - 1.0 * (x[2]+ B * x[1]^2 - 100*B)
  -c(g1,g2)
}

logf <- function(x) {
  -(x[1]^2)/200 - 0.5*(x[2] + B*x[1]^2 - 100*B)^2 + 100 # keeps the shape
}
logf_grad <- function(x) {
  g1 <- -(x[1])/100- 1.0 * (2* B * x[1]) * (x[2]+ B * x[1]^2 - 100*B)
  g2 <- - 1.0 * (x[2]+ B * x[1]^2 - 100*B)
  c(g1,g2)
}

b.M <- matrix(data = c(1, 0, 0, 1), ncol = 2)


# with Metropolis-Hastings

b.M <- matrix(data = c(2, 0.5, 0.5, 1), ncol = 2)

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- MH(logf, b.M, m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}


plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)


# with Hamiltonian MC

# set.seed(4)
samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- HMC(minus_logf, minus_logf_grad, M = b.M, step.size = 0.4, L = 50, q.0 = c(0, 5), m = m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}


plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)


# grid for best params
best.ess <- 0
best.step <- 0
best.L <- 0
best.Q <- 0

for(i in 1:4) {
  if(i == 1) {
    Q <- matrix(data = c(1, 0, 0, 1), ncol = 2)
  } else if(i == 2) {
    Q <- matrix(data = c(2, 0, 0, 1), ncol = 2)
  } else if(i == 3) {
    Q <- matrix(data = c(1, 0.5, 0.5, 1), ncol = 2)
  } else if(i == 4) {
    Q <- matrix(data = c(2, 0.5, 0.5, 1), ncol = 2)
  }
  
  for(step.size in seq(0.2, 1, 0.2)) {
    for(L in seq(10, 50, 10)) {
      cat(step.size, L, "\n")
      # sample
      
      set.seed(67841)
      samples <- NULL
      times <- rep(0, chains)
      for(c in 1:chains) {
        start <- as.numeric(Sys.time())
        sample <- HMC(minus_logf, minus_logf_grad, M = Q, step.size = step.size, L = L, q.0 = c(0, 5), m = m)
        end <- as.numeric(Sys.time())
        
        samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
        times[c] <- end-start
      }
      
      ESS <- array(dim = c(chains, 2))
      ESS.s <- array(dim = c(chains, 2))
      for(c in 1:chains) {
        chain.n <- samples[samples$chain == as.character(c),][,1:2]
        ESS[c,] <- ess(chain.n) # ESS
        ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
      }
      
      if(mean(ESS) > best.ess) {
        best.ess <- mean(ESS)
        best.step <- step.size
        best.L <- L
        best.Q <- Q
      }
    }
  }
}


# with rejection sampling

a <- c(-40, -60)
b <- c(40, 20)
f.max <- 100
M <- f.max*prod(b-a)

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- RS(logf, M, a, b, 2, m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}


plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)


############################################
# c) and d) logistic regression likelihood #
############################################

data <- read.csv("dataset.csv")
idx <-  1:1000# sample(1000, 100)
y <- data$y[idx]
x.2 <- as.matrix(data[idx,1:2])
x.11 <- as.matrix(data[idx,1:11])

loglike <- function(x, w) {
  p <- 1/(1+exp(-x%*%w))
  loss <- y*log(p + 1e-15) + (1-y)*log(1-p + 1e-15)
  return(sum(loss))
}

loglike_grad <- function(x, w) {
  p <- 1/(1+exp(-x%*%w))
  return(c(t(x)%*%(y-p)))
}

c.f.2 <- function(w) {
  return(loglike(x.2, w)+5000)
}
c.f.2.grad <- function(w) {
  return(loglike_grad(x.2, w))
}


x0 <- seq(-1, 3, 0.1)
x1 <- seq(-3, 0, 0.1)
x.graph <- expand.grid(x0,x1)
y.graph <- apply(x.graph, 1, c.f.2)
ggplot()+geom_contour_filled(data=data.frame(x.graph, y=y.graph), aes(Var1, Var2, z=y), alpha = 0.5, colour="black", show.legend=TRUE) # seq(-11000, 0, 100), bins=100


# True

model.2 <- rstanarm::stan_glm(y ~ X1+X2-1, data = data, family = binomial(link = "logit"),
                            algorithm = "sampling", chains = 1, iter = 1000)
samples <- rstan::extract(model.2$stanfit)
w.2 <- array(dim=2)
for(i in 1:2) {
  cat(mean(samples$beta[,i]), "\n")
  w.2[i] <- mean(samples$beta[,i])
}


model.11 <- rstanarm::stan_glm(y ~ .-1, data = data, family = binomial(link = "logit"),
                              algorithm = "sampling", chains = 1, iter = 1000)
samples <- rstan::extract(model.11$stanfit)
w.11 <- array(dim=11)
for(i in 1:11) {
  cat(mean(samples$beta[,i]), "\n")
  w.11[i] <- mean(samples$beta[,i])
}



# with Metropolis-Hastings

Q <- matrix(data = c(0.05, 0, 0, 0.05), ncol = 2)

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- MH(c.f.2, Q, m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}


plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)


# with Hamiltonian MC

Q <- matrix(data = c(0.2, 0, 0, 0.2), ncol = 2)

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- HMC(c.f.2, c.f.2.grad, M = Q, step.size = 0.001, L = 1, q.0 = c(0, 0), m = m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}

plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)


# with rejection sampling

a <- c(-10, -10)
b <- c(10, 10)
f.max <- 1600
M <- f.max*prod(b-a)

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- RS(c.f.2, M, a, b, 2, m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x1=sample[,1], x2=sample[,2], chain=as.character(c)))
  times[c] <- end-start
}


plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 2))
ESS.s <- array(dim = c(chains, 2))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:2]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance") # autocovariance
  if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS[,1])
mean(ESS[,2])
mean(ESS.s[,1])
mean(ESS.s[,2])

mean(samples$x1)
mean(samples$x2)




# d

c.f.11 <- function(w) {
  return(loglike(x.11, w)+5000)
}
c.f.11.grad <- function(w) {
  return(loglike_grad(x.11, w))
}

# with Metropolis-Hastings

Q <- diag(11)/20

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- MH(c.f.11, Q, m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(x=sample, chain=as.character(c)))
  times[c] <- end-start
}
colnames(samples) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "chain")
colMeans(samples[,-12])
w.11


#plot.samples(samples)
plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 11))
ESS.s <- array(dim = c(chains, 11))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:11]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
  acf(chain.n, lag.max = 100, type="covariance", plot=FALSE) # autocovariance
  #if(c < chains) {readline("Press ENTER for next chain...")}
}
mean(ESS)
mean(ESS.s)

# grid
rmse <- function(a, b) {
  return(sqrt(mean((a-b)^2)))
}

best.i <- 0
best.samples <- NULL
best.r <- 100
for(i in seq(10, 100, 10)) {
  Q <- diag(11)/i
  
  samples <- NULL
  for(c in 1:chains) {
    sample <- MH(c.f.11, Q, m)
    samples <- rbind(samples, data.frame(x=sample, chain=as.character(c)))
  }
  r <- rmse(colMeans(samples[,-12]), w.11)
  if(r < best.r) {
    best.r <- r
    best.i <- i
    best.samples <- samples
  }
}


# Hamiltonian MC

Q <- diag(11)*20

samples <- NULL
times <- rep(0, chains)
for(c in 1:chains) {
  start <- as.numeric(Sys.time())
  sample <- HMC(c.f.11, c.f.11.grad, M = Q, step.size = 0.001, L = 5, q.0 = rep(0, 11), m = m)
  end <- as.numeric(Sys.time())
  
  samples <- rbind(samples, data.frame(sample, chain=as.character(c)))
  times[c] <- end-start
}
colnames(samples) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "chain")
colMeans(samples[,-12])
w.11

plot.samples.x1(samples)
plot.samples.x2(samples)

ESS <- array(dim = c(chains, 11))
ESS.s <- array(dim = c(chains, 11))
for(c in 1:chains) {
  chain.n <- samples[samples$chain == as.character(c),][,1:11]
  
  #mcse(as.matrix(chain.n))
  ESS[c,] <- ess(chain.n) # ESS
  ESS.s[c,] <- ESS[c,]/times[c] # ESS/s
}
mean(ESS)
mean(ESS.s)

acf(chain.n, lag.max = 100, type="covariance", plot=FALSE)
