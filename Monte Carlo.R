# Constant gamma(t)
Log_likelihood <- function(gamma){
  
  ## Functions
  Lambda <- function(t){
    return(gamma * (1 - exp(-w*(Time-t))))
  }
  Alpha <- function(t){
    return(gamma * exp(-w*(Time-t)))
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
  
  # m = 1:m.x
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
Simulate <- function(Time, beta, w, gamma){
  y = 0
  tx = 0
  ty = 0
  t.x = c()
  # simulate the first jump in y
  ty = ty + rexp(1, gamma)
  
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
      ty = ty + rexp(1, gamma)
    }
    # both x and y jump after T, end simulation
    else{tx = Time}
  }
  return(t.x)
}

beta = 1
w = 1
epsilon = 1
Time = 10
gamma = 2

gamma_hat = c()
for (i in 1:10000){
  t.x = Simulate(Time, beta, w, gamma)
  gamma_hat[i] = optimize(Log_likelihood, c(0, 5), maximum=TRUE)$maximum
}
hist(gamma_hat, breaks=20, freq=F, main=expression(paste('Histogram of MLE ', hat(gamma))),
     xlab=expression(hat(gamma)))
bias = mean(gamma_hat) - gamma 
mse = mean((gamma_hat - gamma)^2) 



# Linear gamma(t)
Log_likelihood <- function(gamma){
  
  ## Functions
  Gamma <- function(t){
    return(gamma[1] + gamma[2]*t) 
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
  
  # m = 1:m.x
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

beta = 0.1
w = 0.2
epsilon = 0.2
Time = 10
k = 2
gamma.true = c(10, 2)
Gamma.true <- function(t){
  k = length(gamma.true)
  f = 0
  for (i in 1:k){
    f = f + gamma.true[i] * t^(i-1)
  }
  return(f) 
}
gamma_hat = matrix(nrow=1000, ncol=2)

for (i in 1:10){
  t.x = Simulate(Time, beta, w, gamma.true)
  gamma_hat[i,] = optim(c(10, 2), Log_likelihood, method='L-BFGS-B', control=list(fnscale=-1),
                        lower=c(0, 0), upper=c(20, 4))$par
}
hist(gamma_hat[,1], breaks=20)
bias = mean(gamma_hat[,1]) - gamma.true[1]
mse = mean((gamma_hat[,1] - gamma.true[1])^2)

hist(gamma_hat[,2], breaks=20)
bias = mean(gamma_hat[,2]) - gamma.true[2]
mse = mean((gamma_hat[,2] - gamma.true[2])^2)

curve(Gamma.true, ylim = c(0, 50), from = 0, to = Time, lwd=1, col='red',
      main=expression(paste('Lines of MLE ', hat(gamma)(t))), xlab='t', ylab=expression(hat(gamma)(t)))
for (i in 10*(1:300)){
  curve(gamma_hat[i, 1] + gamma_hat[i, 2]*x, from = 0, to = Time, add=T, lwd=0.01)
}
curve(Gamma.true, add=T, from = 0, to = Time, lwd=1, col='red')


