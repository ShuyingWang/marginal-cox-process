### Functions ###

# Simulate Cox process
Simulate <- function(Time, beta, w, gamma){
  y = 0
  tx = 0
  ty = 0
  t.x = c()
  Gamma <- function(t){
    k = length(gamma)
    f = 0
    for (i in 1:k){
      f = f + gamma[i] * t^(i-1)
    }
    return(f) 
  }
  Int_Gamma <- function(t){
    k = length(gamma)
    f = 0
    for (i in 1:k){
      f = f + gamma[i] * t^i / i
    }
    return(f) 
  }
  
  # simulate the first jump in y
  u = rexp(1, 1)
  Dif <- function(t){
    return(Int_Gamma(t) - u)
  }
  if (Dif(Time) < 0){
    ty = Time
  }else{
    ty = uniroot(Dif, c(0, Time), extendInt = 'upX')$root
  }
  
  while(tx < Time){
    
    # simulate a jump in x
    tx = tx + rexp(1, beta + w*y)
    
    # x jumps before y, record this jump
    if (tx < ty & tx < Time){
      t.x = c(t.x, tx)
    }
    # x jumps after y, discard this jump and simulate a new jump in y
    else if (ty < Time){
      y = y + 1
      tx = ty
      u = u + rexp(1, 1)
      Dif <- function(t){
        return(Int_Gamma(t) - u)
      }
      if (Dif(Time) < 0){
        ty = Time
      }else{
        ty = uniroot(Dif, c(0, Time), extendInt = 'upX')$root
      }
    }
    # both x and y jump after T, end simulation
    else{tx = Time}
  }
  return(t.x)
}

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
m0 = c(0, 0, 0)
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
for (i in 500*(1:20)){
  Gamma.sample <- function(t){
    f = 0
    for (j in 1:k){
      f = f + gamma.sample[i,j] * t^(j-1)
    }
    return(f) 
  }
  curve(Gamma.sample, from = 0, to = Time, add=T, lwd=0.1)
}
curve(Gamma.true, add=T, from = 0, to = Time, lwd=1.5, col='red')


