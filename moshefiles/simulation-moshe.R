
#run the functions below, first
#   u=exp(b0*b1*t) : expected introduction rate per year: the mean for the poisson . M= total number of Med species
M=500
b0=0.2
b1=0.05


sim.discovery.data=sim.record(M, b0, b1) #<====================== GENERATE SIMULATED DATA

plot(n.Nativ~t, sim.discovery.data, col='red', ylim=c(0,M+20))  #the number of known nativess
points(n.Inv~t, sim.discovery.data, col='blue')                 #the number of known invasives
lines(I.tot~t, sim.discovery.data)  #the number of invasive species in the Med (known and unknown)
points(I.tot_max.LL~t, sim.discovery.data, col='green', cex=log(n.sps_t),pch=20) #maximum likelihood estimate of I.tot, given binomial
abline(h=M)

#model fitting the 'actual' (simulation) I.tot <---- need to think of a better function than the log-log
lm.res.sim=lm(log(I.tot+0.5)~log(t), sim.discovery.data)
summary(lm.res.sim)
#model fitting the 'extimated' (max-likelihood) I.tot
lm.res=lm(log(I.tot_max.LL+0.5)~log(t), sim.discovery.data, weights=exp(max.LL)^0.25)
summary(lm.res)

plot(log(I.tot+0.5)~log(t), sim.discovery.data, ylim=c(0,8))
abline(lm.res.sim)
points(log(I.tot_max.LL+0.5)~log(t), sim.discovery.data, col='red')
abline(lm.res, col='red')



#the functions for running the above :

multi.pois=function(x){
  r.num=array(0,length(x))
  for(i in 1:length(x)){
    r.num[i]=rpois(1,x[i])
  }
  return(r.num)
}

max.LL.Inv=function(S, N){
  ps=seq(0.01, 0.99, by = 0.01)
  LL=dbinom(x = S, size = N, prob = ps,  log = T)
  max.LL=max(LL)
  id.max=which(LL==max.LL)
  p_max.LL=ps[id.max]
  return(c(max.LL, p_max.LL))
}

obtain.samp=function(I.tot, M, n.Inv, n.Nativ){
  uk.Inv=I.tot-n.Inv
  uk.Nativ=M-n.Nativ
  sps.list=c(rep(1,uk.Nativ), rep(2,uk.Inv), rep(3,n.Nativ), rep(4,n.Inv))   #1=unknow Nativ, 2=unknown Inv, 3=known Nativ, 4-known Inv
  sps.list=sample(sps.list,length(sps.list), replace=F)  #shuffle the list
  return(sps.list)
}

sim.record=function(M, b0, b1){
  t=1:121
  u=exp(b0+b1*t) # expected introduction rate per year: the mean for the poisson 
  sim.data=cbind.data.frame(t,u)
  sim.data$I.t=multi.pois(sim.data$u) # the number of invasives introduced in year t; sampled from poisson distribution with mean u[t]
  sim.data$I.tot=cumsum(sim.data$I.t) #cumulative number of invasives
  sim.data$tot.sps=M+sim.data$I.tot #cumulative total number of species
  
  #sampling effort - how many spcies (in total) are sampled (range prescribed by a min and max number)
  #  - can either increase monotonicaly with time, or be random 
  monotonic=0 # is sampling effort increasing monotonicaly with time?  1=yes, 0=no
  min.sps.sampled=10 #per year
  max.sps.sampled=30
  samp.effort=ceiling(runif(length(t),min=min.sps.sampled, max=max.sps.sampled)) #sampling effort : number of species sampled
  if(monotonic==1) {
    samp.effort=sort(samp.effort, decreasing = F)  #use if monotonic incease in sampling effort
  }
  
  #not all of the t years were sampled
  n.years.0 = 60 # number of years with no sampling
  id.years.0=sample(1:length(t), n.years.0, replace = F)
  samp.effort[id.years.0]=0
  sim.data$samp.effort=samp.effort
  
  # dataframe including only those years that were sampled
  sim.discovery.data=sim.data[which(sim.data$samp.effort>0),]
  sim.discovery.data$n.Nativ=0  #total number of known invasives before sampling of year t
  sim.discovery.data$n.Inv=0    #total number of known natives before sampling of year t
  sim.discovery.data$n.Nativ_t=0 #number of new invasives  sampled on year t
  sim.discovery.data$n.Inv_t=0   #number of new natives  sampled on year t
  sim.discovery.data$n.sps_t=0  # number of new natives & invasives  sampled on year t
  
  
  #obtain sample  
  for (i in 1:(dim(sim.discovery.data)[1]-1)){
    sps.list=obtain.samp(sim.discovery.data$I.tot[i], M, sim.discovery.data$n.Inv[i], sim.discovery.data$n.Nativ[i])
    samp=sample(sps.list, sim.discovery.data$samp.effort[i], replace=F)  #obtain a sample of the prescribed number of species, from the list
    sim.discovery.data$n.Nativ_t[i]=length(which(samp==1))
    sim.discovery.data$n.Inv_t[i]=length(which(samp==2))
    sim.discovery.data$n.Nativ[i+1]=sum(sim.discovery.data$n.Nativ_t)
    sim.discovery.data$n.Inv[i+1]=sum(sim.discovery.data$n.Inv_t)
    sim.discovery.data$n.sps_t[i]=length(which(samp==1 | samp==2))
  }
  
  #remove years in which, by chance, no new species were recorded
  threshold.n=1
  sim.discovery.data=sim.discovery.data[which(sim.discovery.data$n.sps_t>threshold.n),]
  
  j=length(sim.discovery.data$t)
  sim.discovery.data$t1=0
  sim.discovery.data$t1[1]=1
  sim.discovery.data$t2=sim.discovery.data$t
  sim.discovery.data$t1[2:j]=sim.discovery.data$t2[1:(j-1)]+1
  
  sim.discovery.data$max.LL=0
  sim.discovery.data$p_max.LL=0
  sim.discovery.data$I.tot_max.LL=0
  
  for (i in 1:j){
    f.out=max.LL.Inv(S=sim.discovery.data$n.Inv_t[i], N=sim.discovery.data$n.sps_t[i])
    sim.discovery.data$max.LL[i]=f.out[1]  #the likelihood of the most likely probability of sampling a new Invasive, given n.Inv_t and n.sps_t
    sim.discovery.data$p_max.LL[i]=f.out[2] # the corresponding probability
    
    #calculates total number of invasives (known & unknown), given M and p_max.LL
    sim.discovery.data$I.tot_max.LL[i]=(f.out[2]*(M-sim.discovery.data$n.Nativ[i])/(1-f.out[2]))+sim.discovery.data$n.Inv[i]
  }
  return(sim.discovery.data)
}


