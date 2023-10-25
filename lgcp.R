library(lgcp)
library(spatstat)

# Create a data frame
t <- t.x*10
n <- length(t)
x <- runif(n, 0, 5)
y <- runif(n, 0, 5)
win <- owin(xrange=c(0, 5), yrange=c(0, 5))
data <- cbind(x,y,t)
tlim <- c(0,100)

# Convert into a space-time planar point pattern object
xyt <- stppp(list(data = data, tlim = tlim, window = win))
xyt <- integerise(xyt)

# Estimating the spatial and temporal component
den <- density.ppp(xyt, sigma=1)
sar <- spatialAtRisk(den)
mut <- muEst(xyt)
plot(mut)
theta <- thetaEst(xyt, spatial.intensity=sar, temporal.intensity=mut, sigma=0.8, phi=1)

# Run
gfun <- function(Y){
  return(Y)
}
exceed <- exceedProbs(c(1.5, 2))
tmpdr <- tempdir()
lg <- lgcpPredict(  xyt=xyt,
                    T=99,
                    laglength=99,
                    model.parameters=lgcppars(sigma=0.8, phi=1, theta=1.2),
                    cellwidth=1,
                    spatial.intensity=sar,
                    temporal.intensity=mut,
                    mcmc.control=mcmcpars(mala.length=11000,burnin=1000,
                                          retain=10,
                                          adaptivescheme=andrieuthomsh(inith=1,alpha=0.5,C=1,
                                                                       targetacceptance=0.5)),
                    output.control=setoutput(gridfunction=dump2dir(dirname=tmpdr,forceSave=TRUE,lastonly=F),
                                             gridmeans=MonteCarloAverage("exceed")))
ex <- expectation(obj=lg, fun=exceed)



# MCMC diagnostic
plot(hvals(lg)[1000:11000],type="l",xlab="Iteration",ylab="h")
tr <- extract(lg,x=3,y=3,t=1,s=-1)
plot(tr,type="l",xlab="Iteration",ylab="Y")

# Posterior mean intensity
intensity <- intens(lg)
Mean_intens <- c()
for (i in 1:100){
  Mean_intens[i] <- sum(intensity$grid[[i]])
}
mu <- c()
for (i in 1:100){
  mu = c(mu, mut(i-1))
}

# Integrated intensity
plot(t.x, 1:n, type='l')
Int_intens <- c()
Int_mu <- c()
for (i in 1:100){
  Int_intens[i] <- sum(Mean_intens[1:i])
  Int_mu[i] <- sum(mu[1:i])
}
lines((1:100)/10, Int_intens, col='blue')

# Comparison
k <- dim(gamma.sample)[2]
B <- dim(gamma.sample)[1]
b <- 1000
gamma.mean <- colMeans(gamma.sample[b:B,])

# cubic spline version for real example 1
Gamma.mean <- function(t){
  f = 0
  for (j in 1:4){
    f = f + gamma.mean[j] * (t/10)^(j-1)
  }
  f = f + gamma.mean[5] * ((t-40)/10)^3 * (t>40)
  return(f/10)
}
Gamma.integral <- function(t){
  f = 0
  for (j in 1:4){
    f = f + gamma.mean[j] * (t/10)^(j)/(j)
  }
  f = f + gamma.mean[5] * ((t-40)/10)^4/4 * (t>40)
  return((f))
}

# polynomial version for other real examples
Gamma.mean <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * (t/10)^(j-1)
  }
  return(f/10)
}
Gamma.integral <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * (t/10)^(j)/(j)
  }
  return((f))
}

# polynomial version for simulation examples
Gamma.mean <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * (t/10)^(j)/j
  }
  return(beta + w*f/10)
}
Gamma.integral <- function(t){
  f = 0
  for (j in 1:k){
    f = f + gamma.mean[j] * (t/10)^(j+1)/((j+1)*j)
  }
  return(beta*t/10 + w*f)
}

plot(t, 1:n, type='l', main='Count data overlapped with integrated intensity',
     xlab='t', ylab='count')
lines(Int_intens, col='blue')
curve(Gamma.integral, col='red', add=T)



# Posterior probability of exceeding a threshold at time t0
t0 <- 100
gamma_t0 <- c()
for (i in b:B){
  f = 0
  for (j in 1:k){
    f = f + gamma.sample[i, j] * (t0/10)^(j-1)
  }
  gamma_t0[i-b+1] = f/10
}

mean(ex[[t0]][,,1]) # 1.5
mean(ex[[t0]][,,2]) # 2

mean(gamma_t0 > Gamma.mean(t0)*1.5)
mean(gamma_t0 > Gamma.mean(t0)*2)


