source('Functions.R')


### Parameter and hyper-parameter settings ###

# Time interval
Time = 10


## Constant gamma(t)
k = 1
gamma.true = 2

# hyper-parameters
beta = 1
w = 1
epsilon = 1

# Prior mean and sd
m0 = 0
s0 = 100

# Proposal sd
s = 2 



## Linear gamma(t)
k = 2
gamma.true = c(10, 2)

# hyper-parameters
beta = 0.1
w = 0.2
epsilon = 0.2

# Prior mean and sd
m0 = c(0, 0)
s0 = c(100, 10)

# Proposal sd
s = c(3, 1.5) 



## Quadratic gamma(t)
k = 3
gamma.true = c(30, -6, 0.6)

# hyper-parameters
beta = 1
w = 0.1
epsilon = 0.1

# Prior mean and sd
m0 = c(0, 0, 0)
s0 = c(1000, 100, 10)

# Proposal sd
s = c(3, 1.5, 0.3)



## Cubic gamma(t)
k = 4
gamma.true = c(20, 10, -3, 0.2)

# hyper-parameters
beta = 1
w = 0.15
epsilon = 0.1

# Prior mean and sd
m0 = c(0, 0, 0, 0)
s0 = c(1000, 100, 100, 10)

# Proposal sd
s = c(3, 1.5, 0.3, 0.03) 



### Data simulation ###

# Curve of gamma(t)
Gamma.true <- function(t){
  k = length(gamma.true)
  f = 0
  for (i in 1:k){
    f = f + gamma.true[i] * t^(i-1)
  }
  return(f) 
}

# Simulate data
t.x = Simulate(Time, beta, w, gamma.true)



### MCMC ###
acc = 0
B = 10000
gamma = gamma.true
gamma.sample = matrix(nrow=B, ncol=k)
ll = Log_likelihood(t.x, Time, beta, w, gamma, epsilon)
ll.sample = c()
t.grid = seq(0, Time, 0.001)
for (i in 1:B){
  gamma.new = rnorm(k, gamma, s)
  Gamma.new <- function(t){
    f = 0
    for (i in 1:k){
      f = f + gamma.new[i] * t^(i-1)
    }
    return(f) 
  }
  if(min(Gamma.new(t.grid)) >= 0){
    
    ll.new = Log_likelihood(t.x, Time, beta, w, gamma.new, epsilon)
    log.acc = ll.new + sum(dnorm(gamma.new, m0, s0, log = TRUE)) - ll - sum(dnorm(gamma, m0, s0, log = TRUE))
    
    u = runif(1)
    if (log(u) < log.acc){
      gamma = gamma.new
      ll = ll.new
      acc = acc + 1
    }
  }
  gamma.sample[i,] = gamma
  ll.sample[i] = ll
  print(i)
}
acc/B


# Trace plots and histograms
h = 2
par(mfrow=c(k, 2))
for (i in 1:k){
  plot(gamma.sample[,i], type = 'l')
  hist(gamma.sample[,i], freq=F)
  lines(density(gamma.sample[,i], adjust=h))
}


# Posterior curves
par(mfrow=c(1,1))
curve(Gamma.true, ylim = c(0, 50), from = 0, to = Time, lwd=1, col='red',
      main=expression(paste('Curves of Posterior Samples of ', gamma(t))), xlab='t', ylab=expression(gamma(t)))
for (i in 50*(1:20)){
  curve(gamma.sample[i, 1] + gamma.sample[i, 2]*x + gamma.sample[i, 3]*x^2, 
        from = 0, to = Time, add=T, lwd=0.1)
}
curve(Gamma.true, add=T, from = 0, to = Time, lwd=1.5, col='red')


# Posterior mean
par(mfrow=c(1,1))
b <- 1000
gamma.mean <- colMeans(gamma.sample[b:B,])
Gamma.mean <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * (t)^(j-1)
  }
  return((f))
}
Gamma.integral <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * (t)^(j+1)/((j+1)*(j))
  }
  return((beta*t + w*f))
}
curve(Gamma.true, from = 0, to = Time)
curve(Gamma.mean, add=T, col='red')

