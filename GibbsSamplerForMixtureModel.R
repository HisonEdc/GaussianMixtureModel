library(gtools)##must download package from server
n <- 100 ##number of data points
k <- 3 ##number of components
mu <- c(1:k)
lambda <- c(1:k)
pi <- c(1:k)
z <- c(1:k)
beta <- 2
x <- c(1:n)
for(i in 1:k)
{
  mu[i] <- k*i
  lambda[i] <- k*i
  pi[i] <- 1/k }

## simulate data from model
for(i in 1:n)
{
  z[i] <- sample(c(1:k),1,replace = TRUE, prob = pi) ## sample an z according to pi
  x[i] <- rnorm(1,mu[z[i]], 1/sqrt(lambda[z[i]]))
}
hist(x,probability=T,nclass=100,main=" ") ##plot data
N <- 10000 ## number of MCMC iterations
m <- 5## frequency of storing samples
s <- N/m ## number of samples to be stored
mus <- matrix(0,s,k) ## store mu
lambdas <- matrix(0,s,k) ## store lambda
pis <- matrix(0,s,k) ## store pi
betas <- c(1:s) ## store beta
mus[1,] <- mu
lambdas[1,] <- lambda
pis[1,] <- pi
betas[1] <- beta
probs <- c(1:k)

## prior parameters
rg <- (max(x)-min(x))
xi <- rg/2
kappa <- 2/(rg*rg)
alpha <- 2
g <- 0.2
h <- 10/(rg*rg)
delta <- 1
nct <- c(1:k)
ctr <- 2 ## for storing
##main MCMC loop
for(i in 1:N)
{
  ##sample z
  nct <- rep(0,k)
  for(j in 1:n)##loop over data
  {
    probs <- rep(0,k)
    sump <- 0
    for(l in 1:k)
    {
      probs[l] <- pi[l]*sqrt(lambda[l])*exp(-lambda[l]*0.5*(x[j]-mu[l])*(x[j]-mu[l]))
      sump <- sump + probs[l]
    }
    for(l in 1:k)
    {
      probs[l] <- probs[l]/sump
    }
    z[j] <- sample(c(1:k),1,replace = TRUE, prob = probs)
    nct[z[j]] <- nct[z[j]] + 1
  }
  ##sample mu
  for(j in 1:k)
  {
    me <- 0
    for(l in 1:n)
    {
      if(z[l]==j)
      {
        me <- me + x[l]
      }
    }
    me <- (lambda[j]*me + kappa*xi)/(lambda[j]*nct[j] + kappa)
    vari <- 1/(lambda[j]*nct[j] + kappa)
    mu[j] <- rnorm(1,me,sqrt(vari))
  }
  ##sample lambda
  for(j in 1:k)
  {
    me <- 0
    for(l in 1:n)
    {
      if(z[l]==j)
      {
        me <- me + (x[l]-mu[j])*(x[l]-mu[j])
      }
    }
    me <- me*0.5
    me <- me + beta
    al <- alpha + 0.5*nct[j]
    lambda[j] <- rgamma(1,al,me)
  }
  ## update beta
  sump <- sum(lambda)
  beta <- rgamma(1,g+k*alpha,h+sump)
  ## update pi
  for(j in 1:k)
  {
    probs[j] <- delta + nct[j]
  }
  pi <- rdirichlet(1,probs)
  l <- i %% m
  if(l==0)
  {
    mus[ctr,] <- mu
    lambdas[ctr,] <- lambda
    pis[ctr,] <- pi
    betas[ctr] <- beta
    ctr <- ctr + 1
  }
}

##plot MCMC samples
plot(mus[,1],type="l",ylim=c(2,10),ylab="mu",xlab="iteration (every 5th)")
lines(mus[,2],col=2)
lines(mus[,3],col=3)
par(mfrow=c(k,1))
for(i in 1:k)
{
  plot(mus[,i],type="l")
}
par(mfrow=c(k,1))
for(i in 1:k)
{
  plot(lambdas[,i],type="l")
}
par(mfrow=c(k,1))
for(i in 1:k)
{
  plot(pis[,i],type="l")
}
plot(betas,type="l")
