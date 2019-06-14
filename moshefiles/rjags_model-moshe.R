


setwd("C:/Users/mkiflawi/Dropbox/Ideas/Belmaker/")

library(rjags)

modelString="#open quote
model {

for ( i in 1:N ) {
logit_P[i]=a+b*t[i]
P[i]=1/(1+exp(-logit_P[i]))  # <---------------logistic
#P[i]=exp(-a*exp(-b*t[i]))    # <---------------Gompertz
#P[i]=1-exp(-(a*t[i]^b))      # <---------------Weibull

dI[i] ~ dbin(  P[i], dsps[i] )  # p first, then n

res[i] = dI[i] - P[i]*dsps[i] # residual, pr*dsps is the mean
dI.new[i] ~ dbin( P[i], dsps[i] ) 
res.new[i] = dI.new[i] - P[i]*dsps[i]
}

fit <- sum(res[])
fit.new <- sum(res.new[])
a ~ dunif(-10,10)
b ~ dnorm(0, 1) #dunif(-1, 1)

}
" # close quote for modelString
writeLines( modelString , con="lessep_model_pr.txt" ) # write to file

jags.params <- c("a", "b", "P", 'fit','fit.new')  

#initial values for logistic
jags.inits <- function(){
  list("a"=runif(1, -10, 10), "b"=rnorm(1, 0, 1))
}

#initial values for Gompertz
jags.inits <- function(){
  list("a"=runif(1, 1, 10), "b"=runif(1, 0.01, 0.1))
}

#initial values for weibull  -  need to play with these values
jags.inits <- function(){
  list("a"=runif(1, 0, 0.1), "b"=runif(1, 0.1, 10))
}



###sim.discovery.data=sim.record(M, b0, b1) 

jags.data=list(
  t=sim.discovery.data$t,
  N=dim(sim.discovery.data)[1],  #number of observations
  dI=as.integer(sim.discovery.data$n.Inv_t),  #  Invasives added on t
  dsps=as.integer(sim.discovery.data$n.sps_t) # Overall new species sampled on year t
)

#run using rjags
jagsModel = jags.model( file="lessep_model_pr.txt" , data=jags.data , inits=jags.inits,  n.chains=3)  # may need to run several times to initiate
update(jagsModel, n.iter=5000)
samples <- coda.samples(jagsModel, 
                        variable.names=c("a","b", "P", 'fit','fit.new'),
                        n.iter=5000)
summary(samples)
dic.samples(jagsModel,  n.iter=1000)

#compare the maximum-likelihood binomial probabilities and those estimated in rjags
sim.discovery.data$P=as.numeric(summary(samples)$quantiles[1:dim(sim.discovery.data)[1],3])  #the median (50th percentile) of P
plot(P~t, sim.discovery.data, ylim=c(0,1))
plot(P~p_max.LL, sim.discovery.data, ylim=c(0,1))
abline(0,1)

sim.discovery.data$exp_I.tot=sim.discovery.data$n.Inv+(sim.discovery.data$P/(1-sim.discovery.data$P))*(M-sim.discovery.data$n.Nativ)
plot(n.Nativ~t, sim.discovery.data, col='red', ylim=c(0,M+20))
points(n.Inv~t, sim.discovery.data, col='blue')
lines(I.tot~t, sim.discovery.data)
lines(exp_I.tot~t, sim.discovery.data, col='green')
abline(h=M)

#model fitting the 'actual' (simulation) I.tot <---- need to think of a better function than the log-log
lm.res.sim=lm(log(I.tot+0.5)~log(t), sim.discovery.data)
summary(lm.res.sim)

#model fitting the 'extimated' (max. likelihood) I.tot
lm.res=lm(log(exp_I.tot+0.5)~log(t), sim.discovery.data, weights=t^0.5)  #<----note the weights
summary(lm.res)

plot(log(I.tot+0.5)~log(t), sim.discovery.data, ylim=c(0,8))
abline(lm.res.sim)
points(log(exp_I.tot+0.5)~log(t), sim.discovery.data, col='red')
abline(lm.res, col='red')



library(MCMCvis)
MCMCsummary(samples, params = 'a', n.eff = TRUE, digits = 2)  #2.5% to 97.5% gives the 95% credible interval
MCMCsummary(samples, params = 'b', n.eff = TRUE, digits = 2)

MCMCtrace(samples, 
          params = c('a', 'b'),
          ind = T,   #individual chains
          ISB = FALSE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

MCMCplot(samples, params = 'P')

MCMCsummary(samples, round = 2)

