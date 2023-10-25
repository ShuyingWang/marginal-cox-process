### Read data ###
library(reader)
library(data.table)

Timestamp_to_hours <- function(ts){
  ts.h = as.integer(substr(ts, 1, 2))
  ts.min = as.integer(substr(ts, 3, 4))
  ts.sec = as.integer(substr(ts, 5, 6))
  ts.nano = as.integer(substr(ts, 7, 15))
  ts = ts.h + ts.min/60 + ts.sec/3600 + ts.nano/(1e9*3600)
  return(ts)
}

n.readLines('ALL_TRADES_20230403', header=F, skip=0, n = 5)

df = fread(file_path, select=c('Time', 'Exchange', 'Symbol'), 
           colClasses=c('Time'='character'))

IBM = df[df$Symbol=='IBM', c('Time', 'Exchange')]
IBM$Time = Timestamp_to_hours(IBM$Time)

H.IBM = IBM[IBM$Exchange=='H', 'Time']
H.IBM = H.IBM[!is.na(H.IBM)]
H.IBM = (H.IBM-9.5) / (16-9.5)



### Transform data ###
x = H.IBM
Time = 10
x = x*Time
m.x = length(x)
w = 0.2

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

x.t = cbind(sapply(seq(0, Time, 0.01), x_t), seq(0, Time, 0.01))
X.t = cbind(sapply(seq(0, Time, 0.01), X_t), seq(0, Time, 0.01))
plot(X.t[,2], X.t[,1], type='l')

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



### Fit gamma(t) ###

## MCMC
source('Functions.R')
library(truncnorm)
beta = 0
epsilon = 0.1 # for numerical stability when calculating the likelihood
k = 3

set.seed(1)
acc = 0
s = c(3, 1, 0.3) # Need to be adjusted for different dataset
B = 10000
m0 = c(0, 0, 0)
sd0 = c(1000, 100, 100)
gamma = c(10, 0, 0) # An arbitrary initial value
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
  print(i)
}
acc/B



# Trace plots and histograms
h = 5
par(mfrow=c(k, 2))
for (i in 1:k){
  plot(gamma.sample[,i], type = 'l')
  hist(gamma.sample[,i], freq=F)
  lines(density(gamma.sample[,i], adjust=h))
}



# Posterior curves
b = 1000
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
plot(x.t[,2], x.t[,1], xlab='t', ylab='count', type='l', main=expression(paste('Trade counts overlapped with integrated intensity')))
Gamma.integral <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * t^(j)/j
  }
  return(f) 
}
curve(Gamma.integral, from = 0, to = Time, add=T, col='red')


