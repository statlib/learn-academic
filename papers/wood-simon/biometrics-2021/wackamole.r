## simulation based demonstration that you can't bootstrap backwards to get incidence...

england <- c(1,2,0,2,2,0,5,4,1,10,14,20,23,28,40,46,65,63,105,103,149,159,205,264,325,351,359,438,496,574,645,646,697,777,743,726,812,899,790,739,779,718,698,648,685,639,609,570,522,565,484,501,451,437,385,380,343,341,323,312,306,268,251,259,251,266,255,213,202,195,166,183,162,178,170,167,137,155,143,153,149,121,128,115,133,138,120,124,116,91,83,94,108,110,83,86,83,80,73,67,76,49,52,42,58,54,57,48,48,38,44,33,41,50,52,42,33,26)
## ... day 1 March 2nd - lockdown day 23

## combined fatal duration model parameters...
ei2d <- 3.16
si2d <-.42

reps <- 100
n <- sum(england) 
m <- length(england)
d <- rlnorm(n*reps,ei2d,si2d) ## draws from the fatal duration distribution

## create a separate day of death entry for each death, and repeat rep times
day <- rep(rep(1:m-12,england),reps) ## day 1 March 13th - lockdown day 11
lday <- round(day-d) ## `imputed' day of infection

iday <- table(lday) ## supposed incidence


## now go forward...
df <- rlnorm(n*reps,ei2d,si2d) ## draws from the fatal duration distribution
fday <- round(lday + df) ## supposed day of death
dday <- table(fday) ## implied expected deaths

ps <- FALSE
if (ps) postscript("wackamole.eps",width=7,height=5)

par(mar=c(5,5,2,1))
plot(1:m-12,england,xlim=c(-50,120),xlab="day",ylab="deaths/infections",
cex.lab=1.3)
month.axis(start=20,stop=190,origin=73)
ii <- 1:n
for (i in 1:reps) {
  id <- table(lday[ii])
  ii <- ii + n
  lines(as.integer(names(id)),id,col="grey")
}

lines(as.integer(names(iday)),iday/reps,col=1,lwd=2)
lines(as.integer(names(dday)),dday/reps,col=4,lwd=2)

abline(v=11,col=2,lwd=2)
if (ps) dev.off()