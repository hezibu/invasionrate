
library(isotone)

iso.data=sim.discovery.data[which(sim.discovery.data$I.tot_max.LL<400),]
iso.data$wts=exp(iso.data$max.LL)^0.5
iso.reg=gpava(iso.data$t, iso.data$I.tot_max.LL, weights = iso.data$wts, solver = weighted.mean,  ties = "primary", p = NA)  #weighted.fractile weighted.median weighted.mean
#iso.reg=gpava(iso.data$t, iso.data$I.tot_max.LL, weights = iso.data$wts, solver = weighted.fractile,  ties = "primary", p = 0.7) #p=0.5 --> median

plot(n.Nativ~t, iso.data, col='red', ylim=c(0,M+20))
points(I.tot_max.LL~t, iso.data,  cex=-(max.LL),pch=20, col='green')
points(iso.reg$x~iso.reg$z)
points(n.Inv~t, iso.data, pch=3, col='blue')
lines(I.tot~t, iso.data)  #I.tot is the actual (simulated) number of Invasives

plot(iso.reg$x~iso.reg$z)  #$x =iso.reg fitted values; $z=time; $y='observed' values [max.L estimates of I.tot]
points(n.Inv~t, iso.data, pch=3, col='blue')
