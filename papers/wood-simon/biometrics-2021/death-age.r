## https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/weeklyprovisionalfiguresondeathsregisteredinenglandandwales

### SET TO DATA LOICATION
setwd("foo/bar")

## read data 2020 + week 1 of 2021 
dat <- read.table("death-age-2020.txt")
dm <- as.matrix(dat)[,-(1:10)] ## from week 11 13 March
md <- data.frame(d=as.numeric(dm),week=factor(rep(11:53,each=20)),
                age=factor(rep(1:20,43)),a=rep(1:20,43),w=rep(11:53,each=20))
## care home deaths e&w from week 11 13 March 2020
chd <- c(0,2,20,195,826,2050,2794,2423,1666,1660,1090,705,564,369,249,191,169,95,91,69,45,35,40,43,23,17,27,31,38,46,63,106,153,168,280,425,467,571,544,532,602,530,560) 
ad <- colSums(dm[18:20,])
plot(ad);points(chd,col=2)
ph <- 1-chd/ad ## proportion 80+ deaths in hospital
dm1 <- dm
dm1[18:20,] <- t(t(dm1[18:20,])*ph)
md <- data.frame(d=as.numeric(dm1),week=factor(rep(11:53,each=20)),
                age=factor(rep(1:20,43)),a=rep(1:20,43),w=rep(11:53,each=20))


b <- gam(d~s(w,k=30)+s(a,k=8)+ti(w,a,k=c(30,8)),family=nb,data=md)

ps <- FALSE
if (ps) postscript("age-shift.eps",width=11,height=11/3)
par(mfrow=c(1,3),mar=c(5,5,1,1))
c1 <- 1.4;c2=1.3
plot(b,select=1,scheme=1,cex.lab=c1,cex.axis=c2)
plot(b,select=2,scheme=1,cex.lab=c1,cex.axis=c2)
plot(b,scheme=1,select=3,ticktype="detailed",zlim=c(-8,4),theta=135,cex.lab=c1,cex.axis=c2)
if (ps) dev.off()

b0 <- gam(d~ s(w,k=30)+s(a,k=8),family=nb,data=md)


