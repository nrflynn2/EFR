data = colBoth$count
xmins = unique(data) # search over all unique values of data
dat = numeric(length(xmins))
z = sort(data)

for (i in 1:length(xmins)){
  xmin = xmins[i] # choose next xmin candidate
  z1 = z[z>=xmin] # truncate data below this xmin value
  n = length(z1) 
  a = 1+ n*(sum(log(z1/xmin)))^-1 # estimate alpha using direct MLE
  cx = (n:1)/n # construct the empirical CDF
  cf = (z1/xmin)^(-a+1) # construct the fitted theoretical CDF
  dat[i] = max(abs(cf-cx)) # compute the KS statistic
}

D = min(dat[dat>0],na.rm=TRUE) # find smallest D value
xmin = xmins[which(dat==D)] # find corresponding xmin value
z = data[data>=xmin] 
z = sort(z)
n = length(z)
alpha = 1 + n*(sum(log(z/xmin)))^-1 # get corresponding alpha estimate

library(gsl)
library(numDeriv)

# the following code, up to the rpower law, came from this website: http://www.rickwash.com/papers/cscw08-appendix/powerlaw.R

dpowerlaw <-function(x, alpha=2, xmin=1, log=F) {
    if (log)
      log(alpha-1) - log(xmin) - alpha * log(x / xmin)
    else
      ((alpha - 1) / xmin) * ((x / xmin) ^ (-alpha))
}

ppowerlaw <- function(q, alpha=2, xmin=1, lower.tail=T, log.p = F) {
    p <- (q / xmin) ^ (-alpha + 1)
    if (lower.tail)
      p <- 1-p
    if (log.p)
      p <- log(p)
    p
  }

qpowerlaw <- function(p, alpha=2, xmin=1, lower.tail=T, log.p = F) {
    if (!lower.tail)
      p <- 1-p
    if (log.p)
      p <- exp(p)
    xmin * ((1-p) ^ (-1 / (alpha - 1)))
}

rpowerlaw <- function(n, alpha=2, xmin=1) {
    qpowerlaw(runif(n, 0, 1), alpha, xmin)
}

testresult = numeric(2500)
for (i in 1:2500){
  power = rpowerlaw(length(z),alpha,xmin) #randomly generate power law data using the parameters we found
  w = ks.test(z,power) #using KS test to see how good the fit is
  if (w$p.value > 0.10){
    testresult[i] = 1}
  if (w$p.value <= 0.10){
    testresult[i] = 0}
}
sum(testresult)

#FITTING AN EXPONENTIAL:
dat2 = numeric(length(xmins))
z = sort(data)

for (i in 1:length(xmins)){
  xmin = xmins[i] # choose next xmin candidate
  z2 = z[z>=xmin] # truncate data below this xmin value
  n = length(z2) 
  lambda = 1/(mean(z2)- xmin) # estimate lambda using direct MLE
  cx = (1:n)/n # construct the empirical CDF
  cf = 1 - exp(lambda*(xmin - z2)) # construct the fitted theoretical CDF
  dat2[i] = max(abs(cf - cx)) # compute the KS statistic
}

D = min(dat2[dat2>0],na.rm=TRUE) # find smallest D value 
xmin = xmins[which(dat2==D)] # find corresponding xmin value
z = data[data>=xmin] 
z = sort(z)
n = length(z)
lambda = 1/(mean(z) - xmin)

testresult2 = numeric(2500)
for (i in 1:2500){
  expfit = rexp(length(z),lambda) #randomly generate exponential data using the parameters we found
  w1 = ks.test(expfit,z) #using KS test to see how good the fit is
  if (w1$p.value > 0.10){
    testresult2[i] = 1}
  if (w1$p.value <= 0.10){
    testresult2[i] = 0}
}
sum(testresult2)

#REGULAR EXPONENTIAL TEST:
lambda2 = 1/mean(data) 
testresult3 = numeric(length(data))
for (i in 1:2500){
  expfit = rexp(length(data),lambda2) #randomly generate exponential data using the parameters we found
  w2 = ks.test(expfit,data) #using KS test to see how good the fit is
  if (w2$p.value > 0.10){
    testresult3[i] = 1}
  if (w2$p.value <= 0.10){
    testresult3[i] = 0}
}
sum(testresult3)

#LOG NORMAL TEST W/ Xmin:
dat3 = numeric(length(xmins))
z = sort(data)
for (i in 1:length(xmins)){
  xmin = xmins[i] # choose next xmin candidate
  z3 = z[z>=xmin] 
  # truncate data below this xmin value
  n = length(z3) 
  mu = sum(log(z3))/length(z3)
  sigmasq = sum((log(z3) - mu)^2)/length(z3) # estimate lamda using direct MLE
  cx = (1:n)/n # construct the empirical CDF
  cf = pnorm((log(z3) - mu)/sqrt(sigmasq)) # construct the fitted theoretical CDF
  dat3[i] = max(abs(cf - cx)) # compute the KS statistic
}

D = min(dat3[dat3>0],na.rm=TRUE) # find smallest D value
xmin = xmins[which(dat3==D)] # find corresponding xmin value
z = data[data>=xmin] 
z = sort(z)
n = length(z)
mu = sum(log(z))/length(z)
sigmasq = sum((log(z) - mu)^2)/length(z)
testresult4 = numeric(2500)
for (i in 1:2500){
  lognfit = rlnorm(length(data),mean=mu,sd=sqrt(sigmasq)) #randomly generate exponential data using the parameters we found
  w3 = ks.test(lognfit,data) #using KS test to see how good the fit is
  if (w3$p.value > 0.10){
    testresult4[i] = 1}
  if (w3$p.value <= 0.10){
    testresult4[i] = 0}
}
sum(testresult4)

#REGULAR LOG NORMAL TEST:
mu2 = sum(log(data[data>0]))/length(data[data>0])
sigmasq2 = sum((log(data[data>0]) - mu2)^2)/length(data[data>0])
testresult5 = numeric(2500)
for (i in 1:2500){
  lognfit = rlnorm(length(data),mean=mu2,sd=sqrt(sigmasq2)) #randomly generate exponential data using the parameters we found
  w3 = ks.test(lognfit,data) #using KS test to see how good the fit is
  if (w3$p.value > 0.10){
    testresult5[i] = 1}
  if (w3$p.value <= 0.10
  ){
    testresult5[i] = 0}
}
sum(testresult5)