## investigate PCR testing inference reliability...
## Model dP/dt = f - d P where f = exp(Xt b)
## so dP_bj/dt = f X_tj - d P_bj where P_bj is deriv of P wrt b_j
## Y ~ binom(n,P)
## Obviosuly we could write this in terms of numbers in a population of size N,
## but then everything just gets divided by N, so might as well scale it out.


RK4 <- function(g,y0,arg,n,h=1) {
## exectute n full RK4 steps from y0 for dy/dt = g(y,t,arg)
  h2 <- h/2;h6 <- h/6
  y <- y0; t <- 0
  yout <- matrix(y,length(y0),n+1)
  for (i in 1:n) {
    yout[,i] <- y
    k1 <- g(y,t,arg)
    yp <- y + k1 * h2
    k2 <- g(yp,t+h2,arg)
    yp <- y + k2 * h2
    k3 <- g(yp,t+h2,arg)
    yp <- y + k3 * h
    k4 <- g(yp,t + h,arg)
    y <- y + h6*(k1+2*(k2+k3)+k4)
    t <- t + h
  }
  yout[,n+1] <- y
  yout
} ## RK4

test <- function(x2) { ## test prevelance
 x2 <- x2/150
 (0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 * 
            x2)^3 * (1 - x2)^10)/3500
}

g0 <- function(y,t,arg) {
## basic simulator...
  test(t) - arg$delta * y
}

g <- function(y,t,arg) {
## basic semi-parametric incidence model and sensitivities
## Model dP/dt = f - d P where f = exp(Xt b)
## so dP_bj/dt = f X_tj - d P_bj where P_bj is deriv of P wrt b_j
  X <- arg$X; b <- arg$b
  delta <- arg$delta
  i <- round(t/arg$h2)+1
  f <- exp(sum(X[i,]*b)) ## daily prop infected
  dy <- y
  dy[1] <- f - delta * y[1]         ## state
  dy[-1] <- f*X[i,] - delta * y[-1] ## sensitivities
  dy
}


ll <- function(beta,y,arg) {
## penalized log likelihood of incidence from PCR model...
 np <- ncol(arg$X)
 arg$b <- beta
 Sb <- arg$lambda*drop(arg$S %*% beta)
 yo <- RK4(g,rep(0,np+1),arg,arg$n)[,-1]
 p <- yo[1,]
 l <- dbinom(y,arg$N,p,log=TRUE)
 dlp <- y/p - (arg$N-y)/(1-p)
 dl <- drop(yo[-1,] %*% dlp)
 list(ll = sum(l)-sum(beta*Sb)/2, dl = dl - Sb)
}

nl <- function(beta,y,arg) {
  -ll(beta,y,arg)$ll
}

gl <- function(beta,y,arg) {
 -ll(beta,y,arg)$dl
}


ps <- FALSE
if (ps) postscript("PCR-test.eps",width=9,height=3)
par(mfrow=c(1,3),mar=c(5,5,1,1))
set.seed(65)
for (rep in 1:3) {
n <- 100; N <- 400
t <- 1:100
mu <- drop(RK4(g0,0,list(delta=1/10),n))[-1]
y <- rbinom(n,N,mu)

sm <- smoothCon(s(t,k=15),data=data.frame(t=seq(0,n,by=.5)))[[1]]

arg <- list()
arg$X <- sm$X;p <- ncol(arg$X);arg$S <- sm$S[[1]];arg$delta <- .1;
y0 <- arg$delta*y/N;X <- arg$X[1:100*2,]
beta <- as.numeric(coef(lm(log(y0+.01)~X-1)))
arg$b <- beta;  arg$n <- n
arg$N <- N;arg$h2 <- .5;arg$lambda <- 1

for (i in 1:20) {
  er <- optim(beta,nl,gr=gl,method="BFGS",hessian=TRUE,y=y,arg=arg)
  beta <- er$par
  f.hat <- exp(arg$X%*%beta)
  arg$lambda <- (13/arg$lambda - sum(diag(solve(er$hessian,arg$S))))/(beta%*%arg$S%*%beta)*arg$lambda
}
V <- solve(er$hessian)
se <- rowSums(arg$X*(arg$X%*%V))^.5
tt <- 1:length(f.hat)*.5
plot(tt,f.hat,type="l",ylim=c(0,0.004),ylab="incidence",xlab="day")
points(.1*y/N,col="grey")
lines(tt,exp(log(f.hat)+2*se),lty=2)
lines(tt,exp(log(f.hat)-2*se),lty=2)
lines(t,test(t),col=2)
}
if (ps) dev.off()

