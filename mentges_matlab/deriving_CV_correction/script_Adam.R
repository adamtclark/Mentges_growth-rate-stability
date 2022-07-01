## Adams script for CV correction

setwd("/Users/am41xite/Nextcloud/Project_III/CV_correction")
source("functions_Adam.R")

rlst<-seq(1,2, length=20)
varout<-numeric(length(rlst))
cvout<-numeric(length(rlst))

K<-10

for(i in 1:length(rlst)) {
  #               f:1/lambda  d:mu d_sd: sigma sf:stepwidth
  d<-symdyn(r=rlst[i], f=2, d=0, d_sd=0.1, sf=0.1, tmax=100, stochd = TRUE, stocht = TRUE, as.matrix = TRUE)
  x<-d[d[,"time"]>20,"state"]
  t<-d[d[,"time"]>20,"time"]
  N<-x+K
  
  varout[i]<-var(N)
  cvout[i]<-sd(N)/mean(N)
}

#plot(t,N, type="l")

plot(rlst , varout)
mu<-0
sigma<-0.1
lambda<-0.5

#var(x) = (mu^2+sigma^2)*lambda/(2*rlst)
lines(rlst, (mu^2+sigma^2)*lambda/(2*rlst))
#So, that seems to work



#now, try for CV
#CV = sd(N)/mean(N)
#var(x) = var(N) = (CV*mean)^2
#and here, mean = K = 10

plot((cvout*K)^2, varout)
abline(a=0, b=1, lty=3)

plot((cvout*K)^2, (mu^2+sigma^2)*lambda/(2*rlst))
abline(a=0, b=1, lty=3)
#good, those match



#so, now we know
# (CV*mean(N))^2 == var(x) == (mu^2+sigma^2)*lambda/(2*rlst)
# CV^2 == (mu^2+sigma^2)*lambda/(2*rlst)/mean(N)^2
# CV^2*rlst == (mu^2+sigma^2)*lambda/(2)/mean(N)^2 == constant
## mu = 0, sigma = 0.2, lambda=1, mean(N)=K
plot(cvout^2*rlst)
abline(h=(0^2+0.1^2)*0.5/(2)/K^2, col="yellow", lwd = 4)
# abline(h=(mu^2+sigma^2)*lambda/(2)/mean(N)^2, col="red")
#seems to work!