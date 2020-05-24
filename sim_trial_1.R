#Creating Solow and Costello simulations:

#poisson introduction rate: N(t) = mu_t * t 
# Y_t  =  sum(s=1 to t) Y_st 
#mu_t = exp(b0 + b1*t)
#pi_st = exp[g0 + g1 * t + g2* ( exp t - s )]/ 1 + [g0 + g1 * t + g2* ( exp t - s )]
#to go from model to sim - get probability from count_p, get n from count_m, then use for bernoulli distribution?

sim.parameters = set_constant_params(c(-1.0,0.00,-1.1,0.00000,0.0000000),
                                 parameters = c("beta0","beta1","gama0","gama1","gama2"))

sim.results = sim(params = sim.parameters,N = 60)

sim.guess <- set_params_to_optimize(c(-1.0,0.00,-0.8,0.00000,0.0000000),
                                 parameters = c("beta0","beta1","gama0","gama1","gama2"))


optimized <- matrix(NA,10,5)

for (i in seq_len(10)){
  print(i)
optimized[i,] <- as.numeric(optimx::optimr(fn = count_log_like,par = sim.guess, first_record_data = sim.results,const = constant)$par)
}

optim_nm(fun = count_log_like,k = 5,start = sim.guess,trace = T)


C1_sim <- count_lambda(params = optimized$par,N = length(sim.results),const = constant)

optimized$par

plot <- data.frame(t(rbind(1:145,sim.results,C1_sim)))
plot <- gather(plot,key = type,value = obs,2:3)
plot <- plot %>% 
  group_by(type) %>% arrange(V1) %>% mutate(cs = cumsum(obs))

plot %>% ungroup() %>%
  #mutate(type = fct_relevel(type,c("num_discov","C2","C1")))%>% 
  ggplot()+
  aes(x = V1, y = cs,group = type,linetype = type)+
  geom_line()+
  theme_classic()+theme(legend.position = "none")
