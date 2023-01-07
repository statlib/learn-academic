
## now use parametric bootstrap for overall distribution uncertainty
nll <- function(theta,x) { ## log normal negative log lik
  -sum(dlnorm(x,theta[1],theta[2],log=TRUE))
} ## nll
set.seed(6)

mu <- 20.2;sd <- 11.6 ## Linton parameters...
lsd1 <- sqrt(log(sd^2/mu^2+1)); lmu1 <- log(mu) - lsd1^2/2
n.rep <- 100
th1 <- th <- matrix(0,n.rep,2)
th0 <- c(lmu1,lsd1) ## initial values
ns <- 10000 ## number of points for approx KL optimization
inc <- exp(qnorm((1:ns-.05)/ns,1.63,.5)) ## 'sample' from infection to onset dist
for (i in 1:n.rep) { ## bootstrap loop
  ## parameteric bootstrap from the 3 literature models...
  xl <- rlnorm(34,lmu1,lsd1)  ## Linton
  xw <- rgamma(41,shape=20^2/100,scale=100/20) ## Wu
  xv <- rgamma(24,shape=17.8/4,scale=4) ## Verity
  ## Fit the combined model to the sample...
  er <- optim(th0,nll,method="BFGS",x=c(xl,xw,xv))
  th[i,] <- er$par
  ## Fit lognormal model to infection to death distribution
  ## by approximate KL minimization...
  xs <- rlnorm(ns,th[i,1],th[i,2])+inc
  er <- optim(th[i,],nll,method="BFGS",x=xs)
  th1[i,] <- er$par
}
thm <- colMeans(th)   ## mean params onset-to-death
thm1 <- colMeans(th1) ## mean params infection-to-death

## plot the results
ps <- FALSE
if (ps) postscript("ddist.eps",height=2.5,width=9)
x0 <- seq(0,70,by=.1);
lmu <- 2.891236; lsd <- 0.5560864 ## CHESS parameters
chess.th <- c(3.1860616,0.4434512)
par(mfrow=c(1,3),mar=c(5,5,1,1))
plot(x0,dlnorm(x0,lmu,lsd),type="l",ylim=c(0,0.06),xlab="day",ylab="probability",col=4)
lines(x0,dlnorm(x0,lmu1,lsd1),lty=2)
lines(x0,dgamma(x0,shape=17.8/4,scale=4),lty=3)
lines(x0,dgamma(x0,shape=20^2/100,scale=100/20),lty=4)

d0 <- dlnorm(x0,thm[1],thm[2]) 
plot(x0,d0,col=2,type="l",ylim=c(0,0.06),ylab="probability",xlab="day")
for (i in 1:n.rep) {
  d <- dlnorm(x0,th[i,1],th[i,2])
  lines(x0,d,col="grey")
}

lines(x0,d0,col=2,lwd=2)
lines(x0,dlnorm(x0,lmu,lsd),lwd=1,col=4)

d0 <- dlnorm(x0,thm1[1],thm1[2]) 
plot(x0,d0,col=2,type="l",ylim=c(0,0.06),ylab="probability",xlab="day")
for (i in 1:n.rep) {
  d <- dlnorm(x0,th1[i,1],th1[i,2])
  lines(x0,d,col="grey")
}
lines(x0,d0,col=2,lwd=2)
lines(x0,dlnorm(x0,chess.th[1],chess.th[2]),col=4,lwd=1)

if (ps) dev.off()

## make sure setwd is used to set working directory to code location
save(th1,file="o2d-params.rda")
#load("o2d-params.rda")
