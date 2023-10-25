### Prepare Dataset ###

## Load in data
library(lgcp)
load('v52i04.RDA')
x <- xyt$x
y <- xyt$y
t <- xyt$t
win <-xyt$window
tlim <- c(0,100)
data <- cbind(x,y,t)



## Cut and pick windows
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
win1 = data[,1]>=xcut[i] & data[,1]<=xcut[i+1] & data[,2]>=ycut[j] & data[,2]<=ycut[j+1]

# Window 2
i = 1
j = 1
win2 = data[,1]>=xcut[i] & data[,1]<=xcut[i+1] & data[,2]>=ycut[j] & data[,2]<=ycut[j+1]

# Window 3
i = 3
j = 3
win3 = data[,1]>=xcut[i] & data[,1]<=xcut[i+1] & data[,2]>=ycut[j] & data[,2]<=ycut[j+1]



## Transform data
x = data[win1, 3]
x = x/10
Time = 10
m.x = length(x)
w = 0.194 # window 1
#w = 0.26 # window 2
#w = 0.197 # window 3

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



### Fit gamma(t) ###

## MCMC
source('Functions.R')
library(truncnorm)
beta = 0
epsilon = 0.1

# Window 1
k = 5
s = c(2.5, 1.5, 0.5, 0.1, 0.007)*0.6

# Window 2
#k = 4
#s = c(2, 1.5, 1, 0.15)*0.3

# window 3
#k= 5
#s = c(3, 1.5, 1, 0.3, 0.01)*0.2 

B = 10000
acc = 0
m0 = c(0, 0, 0, 0, 0)
sd0 = c(1000, 100, 100, 10, 10)
gamma = c(20, 0, 0, 0, 0)
gamma.sample = matrix(nrow=B, ncol=k)
ll = Log_likelihood(t.x, Time, beta, w, gamma, epsilon)
ll.sample = c()
gamma.new = c()
t.grid = seq(0, Time, 0.001)
for (i in 1:B){
  gamma.new[1] = rtruncnorm(1, 0, Inf, gamma[1], s[1])
  gamma.new[2:k] = rnorm(k-1, gamma[2:k], s[2:k])
  Gamma.new <- function(t){
    f = 0
    for (i in 1:k){
      f = f + gamma.new[i] * t^(i-1)
    }
    return(f) 
  }
  if(min(Gamma.new(t.grid)) >= 0){
    
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
b = 1000
par(mfrow=c(2,1))
plot(gamma.sample[b:B,1], type='l')
plot(gamma.sample[b:B,2], type='l')
plot(gamma.sample[b:B,3], type='l')
plot(gamma.sample[b:B,4], type='l')
plot(gamma.sample[b:B,5], type='l')
plot(ll.sample, type='l')


h = 5
par(mfrow=c(2,2))
hist(gamma.sample[b:B,1], breaks=20, freq=F)
lines(density(gamma.sample[b:B,1], adjust=h))
hist(gamma.sample[b:B,2], breaks=20, freq=F)
lines(density(gamma.sample[b:B,2], adjust=h))
hist(gamma.sample[b:B,3], breaks=20, freq=F)
lines(density(gamma.sample[b:B,3], adjust=h))
hist(gamma.sample[b:B,4], breaks=20, freq=F)
lines(density(gamma.sample[b:B,4], adjust=h))
hist(gamma.sample[b:B,5], breaks=20, freq=F)
lines(density(gamma.sample[b:B,5], adjust=h))



# Posterior curves
b = 100
gamma.mean = colMeans(gamma.sample[b:B,])
gamma.mode = c(estimate_mode(gamma.sample[b:B,1], h), estimate_mode(gamma.sample[b:B,2], h),
               estimate_mode(gamma.sample[b:B,3], h), estimate_mode(gamma.sample[b:B,4], h),
               estimate_mode(gamma.sample[b:B,5], h))
par(mfrow=c(1,1))
curve(gamma.mode[1] + gamma.mode[2]*x + gamma.mode[3]*x^2 + gamma.mode[4]*x^3 + gamma.mode[5]*x^4, from = 0, to = Time,
      xlab='t', ylab=expression(hat(gamma)(t)), main=expression(paste('The fitted ', gamma(t), ' curve')))
curve(gamma.mean[1] + gamma.mean[2]*x + gamma.mean[3]*x^2 + gamma.mean[4]*x^3 + gamma.mean[5]*x^4, col='blue', add=T)

par(mfrow=c(1,1))
plot(x.t[,2], x.t[,1], xlab='t', ylab='x', type='l',
     main=expression(paste('The original data with fitted mean function')))
curve(gamma.mean[1]*x + gamma.mean[2]*x^2/2 + gamma.mean[3]*x^3/3 + gamma.mean[4]*x^4/4 + gamma.mean[5]*x^5/5,
      add=T, col='red')


# AIC and BIC
ll = Log_likelihood(t.x, Time, beta, w, gamma.mean, epsilon)
BIC = k*log(length(t.x)) - 2*ll
AIC = k*2 - 2*ll
