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

# Log likelihood with cubic spline gamma(t)
Log_likelihood_spline <- function(t.x, Time, beta, w, gamma, epsilon=1, t_knot){
  
  ## Functions
  Gamma <- function(t){
    k = 4
    f = 0
    for (i in 1:k){
      f = f + gamma[i] * t^(i-1)
    }
    f = f + gamma[5] * (t-t_knot)^3 * (t>t_knot)
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

# Posterior mode estimation
estimate_mode <- function(x, h=1){
  d <- density(x, adjust=h)
  return(d$x[which.max(d$y)])
}


