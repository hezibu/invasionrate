sim.parameters.timespan.3 = set_constant_params(c(-1.1,0.014,-1.46,0,0),
                                              parameters = c("beta0","beta1","gama0","gama1","gama2"))

sim.matrix.3 <- replicate(expr = sim(params = sim.parameters.timespan.3, N = 150),n = 100)

params.to.optimize.3 = set_params_to_optimize(c(-1.1106, 0.0135, -1.4534),
                                                parameters = c("beta0","beta1","gama0"))

constant <- set_constant_params(c(0,0),parameters = c("gama1","gama2"))

# library(pacman)
# p_load(pbapply,optimx,purrr)


get_parameters_estimate <- function(matrix,params,const){ # matrix of [timeseries , simulations]
  library(purrr)
  estimates <- apply(matrix,2,function(x)
    (optim(fn = count_log_like,
                     par = params,
                     first_record_data = x,
                     const = const, hessian = T)))
  return(estimates)
}



get_parameters_specified_timespan <- function(matrix,params,const,timespan){
  take_nested_timespan(matrix,timespan) %>% 
    get_parameters_estimate(params,const) %>% 
    return()
}


lapply(list(50),function(i) get_parameters_specified_timespan(sim.matrix.3[,1:2],
                                                             params = params.to.optimize.3,
                                                             const = sim.parameters.timespan.3,
                                                             timespan = 50))


parameters.results.timespans <- list()
timespan.intervals <- seq(10,150,by=10)
for (i in seq_along(timespan.intervals)){
  print(i)
  parameters.results.timespans[[i]] <- get_parameters_specified_timespan(sim.matrix,
                                                                         sim.parameters.timespan,
                                                                         timespan.intervals[i])
}

parameters.results.timespans.3 <- pblapply(timespan.intervals,
                                         function(i) get_parameters_specified_timespan(sim.matrix.3,params = sim.parameters.timespan.3[1:3],
                                                                                       const = sim.parameters.timespan.3[4:5],
                                                                                       timespan =  i))


did_optim_succeed <- function(opt,params){
  simulated.params <- as.data.frame(params) %>% 
    rownames_to_column()
  colnames(simulated.params) <- c("par","true")
  
  out <- t(get_parameter_se(opt)) %>% 
    as_tibble() %>% 
    mutate(par = names(opt$par)) %>% 
    left_join(simulated.params,by = "par") %>% 
    mutate(correct = (true <= `95% Upper bound` & 
                        true >= `95% Lower bound`)) %>% 
    select(par,`95% Lower bound`,ML,`95% Upper bound`,true,correct)
  return(out)
}

did_optim_succeed(parameters.results.timespans.3[[15]][[100]],sim.parameters.timespan.3)

onefifty.success <- do.call(bind_rows,
        lapply(parameters.results.timespans.3[4], function(x) lapply(x,function(y) did_optim_succeed(y,sim.parameters.timespan.3))) %>% 
  unlist(recursive = F))

onefifty.success %>% group_by(par) %>% nest() %>% 
  mutate(percent = map(.x = data, function(x) sum(x$correct)/nrow(x))) %>% 
  select(par,percent) %>% 
  unnest()

onefifty.success %>% filter(par == "beta1") %>% 
  mutate(index = as.numeric(rownames(.))) %>%
  ggplot()+
  geom_errorbar(aes(ymin = `95% Lower bound`,ymax = `95% Upper bound`, x = index,color = correct))+
  theme_classic()+
  coord_cartesian(ylim = c(-0.05,0.05))+
  geom_hline(yintercept = 0.014)

  