sim.parameters.constant.rate = set_constant_params(c(3,0,-7,0,0),
                                                parameters = c("beta0","beta1","gama0","gama1","gama2"))


sim.matrix.const.rate <- replicate(sim(params = sim.parameters.constant.rate,N = 150),n = 100)

guess.constant.rate <- set_params_to_optimize(numeric_vector = c(0,-5),parameters = c("beta1","gama0"))
guess.constant.rate.all <- set_params_to_optimize(c(-1,0.015,-2,0,0))

parameters.results.constant.rate <- pblapply(150, function(i) get_parameters_specified_timespan(sim.matrix.const.rate,
                                                                                           params = guess.constant.rate.all,
                                                                                           const = sim.parameters.constant.rate,
                                                                                           timespan =  i))

parameters.results.constant.rate <- parameters.results.constant.rate %>% unlist(recursive = F)
                             
results.did.succeed.constant.rate <- 
  lapply(parameters.results.constant.rate, function(x) try(did_optim_succeed(x,sim.parameters.constant.rate),silent = T))

j <- 0
results.did.succeed.constant.rate.clean <- vector(mode = "list")
for (i in seq_len(length(results.did.succeed.constant.rate))){
  if (is.tibble(results.did.succeed.constant.rate[[i]])){
    results.did.succeed.constant.rate.clean[[i]] <- results.did.succeed.constant.rate[[i]]
    
  }
  else {
    j <- j+1 #count the times getting SE did not work
  }
}

results.did.succeed.constant.rate.clean <- do.call(bind_rows,results.did.succeed.constant.rate.clean)

results.did.succeed.constant.rate.clean %>% group_by(par) %>% nest() %>% 
  mutate(percent = map(.x = data, function(x) sum(x$correct,na.rm = T)/(nrow(x)+j))) %>% 
  select(par,percent) %>% 
  unnest()

colors = setNames( c("red","blue","grey"),
                   levels(factor((results.did.succeed.constant.rate.clean$correct))))

plot.data <- results.did.succeed.constant.rate.clean %>% 
  group_by(par) %>%
  nest() %>%
  mutate(plot = map(data, function(x)
    x %>% mutate(index = as.numeric(rownames(x))))) %>% select(par,plot) %>%  unnest() %>% 
  group_by(par) %>% 
  nest() %>% mutate(plot = map2(data,par, function(x,y) {
    ggplot(x)+
      geom_errorbar(data = x, aes(ymin = `95% Lower bound`,ymax = `95% Upper bound`, x = index,color = correct))+
      geom_point(data = subset(x,is.na(correct)), aes(y = ML, x = index),color = "red")+
      theme_classic()+
      coord_cartesian(ylim = c(mean(x$true) - 1, 
                               mean(x$true) + 1))+
      geom_hline(yintercept = mean(x$true))+
      scale_color_manual(values = colors)+
      ggtitle({y})
  })) %>% 
  select(plot)

coefs_plot <- gridExtra::arrangeGrob(grobs = plot.data$plot,ncol = 1)
ggsave(coefs_plot,filename = "coefs_plot.png",
       width = 8.45, height =  20.72)


parameters.results.constant.rate.fixed <- get_parameters_specified_timespan(sim.matrix.const.rate,
                                                                      params = guess.constant.rate,
                                                                      const = sim.parameters.constant.rate,
                                                                      timespan =  150)

#parameters.results.constant.rate.fixed <- parameters.results.constant.rate.fixed %>% unlist(recursive = F)

results.did.succeed.constant.rate.fixed <- 
  lapply(parameters.results.constant.rate.fixed, function(x) try(did_optim_succeed(x,sim.parameters.constant.rate),silent = T))

results.did.succeed.constant.rate.clean.fixed <- do.call(bind_rows,results.did.succeed.constant.rate.fixed)


plot.data.fixed <- results.did.succeed.constant.rate.clean.fixed %>% 
  group_by(par) %>%
  nest() %>%
  mutate(plot = map(data, function(x)
    x %>% mutate(index = as.numeric(rownames(x))))) %>% select(par,plot) %>%  unnest() %>% 
  group_by(par) %>% 
  nest() %>% mutate(plot = map2(data,par, function(x,y) {
    ggplot(x)+
      geom_errorbar(data = x, aes(ymin = `95% Lower bound`,ymax = `95% Upper bound`, x = index,color = correct))+
      geom_point(data = subset(x,is.na(correct)), aes(y = ML, x = index),color = "red")+
      theme_classic()+
      geom_hline(yintercept = mean(x$true))+
      scale_color_manual(values = colors)+
      ggtitle({y})
  })) %>% 
  select(plot)

comparison.data <- data.frame(two_parameters = unlist(lapply(parameters.results.constant.rate.fixed, function(x) x$value)),
                              five_parameters = unlist(lapply(parameters.results.constant.rate, function(x) x$value)))

comparison.data <- comparison.data %>% mutate(delta = five_parameters - two_parameters)
                              
                              