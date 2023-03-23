### Prepare Dataset ###

## Load in data
load('v52i04.RDA')
x <- xyt$x
y <- xyt$y
t <- xyt$t
data <- cbind(x,y,t)



## Cut and pick windows (rectangle areas in the map)
xl = min(data[,1])
xu = max(data[,1])
xr = xu - xl

yl = min(data[,2])
yu = max(data[,2])
yr = yu - yl

nx = 5
ny = 3
xcut = seq(xl, xu, xr/nx)
ycut = seq(yl, yu, yr/ny)

par(mfrow = c(ny, nx))
for (j in ny:1){
  for(i in 1:nx){
    window = data[data[,1]>=xcut[i] & data[,1]<=xcut[i+1] & data[,2]>=ycut[j] & data[,2]<=ycut[j+1], 3]
    plot(window, 1:length(window), type="l", main = paste('x=',i, ' y=', j))
  }
}

# Window 1
i = 5
j = 2
x = data[data[,1]>=xcut[i] & data[,1]<=xcut[i+1] &
                data[,2]>=ycut[j] & data[,2]<=ycut[j+1], 3]

# Window 2
i = 1
j = 1
x = data[data[,1]>=xcut[i] & data[,1]<=xcut[i+1] &
           data[,2]>=ycut[j] & data[,2]<=ycut[j+1], 3]

# Window 3
i = 3
j = 3
x = data[data[,1]>=xcut[i] & data[,1]<=xcut[i+1] &
           data[,2]>=ycut[j] & data[,2]<=ycut[j+1], 3]



## Transform data
x = x/10
Time = 10
m.x = length(x)
#w = 0.194 # window 1
#w = 0.26 # window 2
w = 0.197 # window 3

x_t <- function(t){
  return(sum(x<=t))
}
X_t <- function(t){
  xt = x[x<=t]
  return(w * sum(t - xt))
}
Integrate_X <- function(t){
  xt = x[x<=t]
  return(w * sum((t - xt)^2) / 2)
}

par(mfrow = c(1, 1))
x.t = cbind(sapply(seq(0, Time, 0.01), x_t), seq(0, Time, 0.01))
X.t = cbind(sapply(seq(0, Time, 0.01), X_t), seq(0, Time, 0.01))
plot(x.t[,2], x.t[,1], xlab='t', ylab='x', type='l', main=expression(paste('The original data overlapped with fitted ', gamma(t))))
lines(X.t[,2], X.t[,1], col='red')

m = floor(X_t(Time))
t.reach = c()
t.jump = c()
Jump_time <- function(t){
  return(X_t(t) - 1)
}
t.reach[1] = uniroot(Jump_time, c(0, Time), extendInt = 'upX')$root
t.jump[1] = t.reach[1] - Integrate_X(t.reach[1])
for (i in 2:m){
  Jump_time <- function(t){
    return(X_t(t) - i)
  }
  t.reach[i] = uniroot(Jump_time, c(0, Time), extendInt = 'upX')$root
  area = Integrate_X(t.reach[i]) - Integrate_X(t.reach[i-1]) - (i-1)*(t.reach[i]-t.reach[i-1])
  t.jump[i] = t.reach[i] - area
}

X_t <- function(t){
  return(sum(t.jump<=t))
}
X.t = cbind(sapply(seq(0, Time, 0.001), X_t), seq(0, Time, 0.001))
lines(X.t[,2], X.t[,1], col='red')
t.x = t.jump
length(t.x)
m.x




### Fit gamma(t) ###

## Functions

# Log likelihood with polynomial gamma(t)
Log_likelihood <- function(t.x, Time, beta, w, gamma, epsilon=1){
  
  ## Functions
  Gamma <- function(t){
    k = length(gamma)
    f = 0
    for (i in 1:k){
      f = f + gamma[i] * t^(i-1)
    }
    return(f) 
  }
  Lambda <- function(t){
    return(Gamma(t) * (1 - exp(-w*(Time-t))))
  }
  Alpha <- function(t){
    return(Gamma(t) * exp(-w*(Time-t)))
  }
  
  ## Integrals
  lambda = integrate(Lambda, lower = 0, upper = Time)$value
  
  m.x = length(t.x)
  t.x = sort(t.x, decreasing=TRUE)
  alpha = c()
  for (i in 1:m.x){
    alpha[i] = integrate(Alpha, lower = 0, upper = t.x[i])$value
  }
  
  ## Coefficients
  
  Recursion <- function(i){
    return(C.m_1[i+1] * choose(m-i-1, j-i-1))
  }
  
  # m = 0, j = 0
  C.m = 1
  
  # m = 1:m.x (At each iteration, multiply the coefficients by a small constant epsilon to prevent them from becoming infinity.
  #            The likelihood calculated in this way is proportional to the original likelihood.)
  for (m in 1:m.x){
    
    C.m_1 = c(C.m, 0)
    C.m = rep(0, m+1)
    
    # j = 0
    C.m[1] = C.m_1[1] * epsilon
    
    # j = 1:m
    for (j in 1:m){
      C.m[j+1] = sum(Recursion(0:(j-1)))
    }
    C.m[2:(m+1)] = (C.m[2:(m+1)] * alpha[m] + C.m_1[2:(m+1)]) * epsilon
  }
  
  ## Log likelihood
  p = 0
  for (j in 0:m.x){
    p = p + C.m[j+1] * w^j * beta^(m.x-j) * exp(-beta*Time -lambda)
  }
  l = log(p)
  
  return(l)
}

# Posterior mode estimation
estimate_mode <- function(x, h=1){
  d <- density(x, adjust=h)
  return(d$x[which.max(d$y)])
}



## MCMC
library(truncnorm)
beta = 0
epsilon = 0.1 # for numerical stability when calculating the likelihood
Time = 10
k = 5

acc = 0
s = c(3, 1.5, 0.5, 0.05, 0.005) # Need to be adjusted for different dataset
B = 10000
m0 = c(0, 0, 0, 0, 0)
sd0 = c(1000, 100, 100, 10, 1)
gamma = c(20, 0, 0, 0, 0) # An abtrary initial value
gamma.sample = matrix(nrow=B, ncol=k)
ll = Log_likelihood(t.x, Time, beta, w, gamma, epsilon)
ll.sample = c()
gamma.new = c()
t.grid = seq(0, Time, 0.001)
for (i in 1:B){
  gamma.new[1] = rtruncnorm(1, 0, Inf, gamma[1], s[1]) # The first coefficient must be positive
  gamma.new[2:k] = rnorm(k-1, gamma[2:k], s[2:k])
  Gamma.new <- function(t){
    f = 0
    for (i in 1:k){
      f = f + gamma.new[i] * t^(i-1)
    }
    return(f) 
  }
  if(min(Gamma.new(t.grid)) >= 0){ # Reject immediately if gamma(t) has negative values in [0, T]
    
    ll.new = Log_likelihood(t.x, Time, beta, w, gamma.new, epsilon)
    log.acc = ll.new + sum(dnorm(gamma.new, m0, sd0, log = TRUE)) + log(dtruncnorm(gamma[1], 0, Inf, gamma.new[1], s[1])) - 
      ll - sum(dnorm(gamma, m0, sd0, log = TRUE)) - log(dtruncnorm(gamma.new[1], 0, Inf, gamma[1], s[1]))
    
    u = runif(1)
    if (log(u) < log.acc){
      gamma = gamma.new
      ll = ll.new
      acc = acc + 1
    }
  }
  gamma.sample[i,] = gamma
  ll.sample[i] = ll
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
gamma.mean = colMeans(gamma.sample[b:B,])
gamma.mode = c()
for (i in 1:k){
  gamma.mode = c(gamma.mode, estimate_mode(gamma.sample[b:B,i], h))
}
par(mfrow=c(1,1))
Gamma.mean <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * t^(j-1)
  }
  return(f) 
}
Gamma.mode <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mode[j] * t^(j-1)
  }
  return(f) 
}
curve(Gamma.mean, from = 0, to = Time,
      xlab='t', ylab=expression(hat(gamma)(t)), main=expression(paste('The fitted ', gamma(t), ' curve')))
curve(Gamma.mode, add=T, col='blue')



par(mfrow=c(1,1))
plot(x.t[,2], x.t[,1], xlab='t', ylab='x', type='l', main=expression(paste('The original data overlapped with fitted ', gamma(t))))
Gamma.integral <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * t^(j)/j
  }
  return(f) 
}
curve(Gamma.integral, from = 0, to = Time, add=T, col='red')


