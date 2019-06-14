heatmap.data.list <- vector(mode = "list",length = length(seq(-0.01,0.03,by = 0.0025))*length(seq(-3,3,by = 0.5)))
#this will store all the data

i <- 1
for (beta1 in seq(-0.01,0.03,by = 0.0025)){
  for (gama0 in seq(-3,3,by = 0.5)){
    simulated <- lapply(1:100, function(i) {
      simulated <- sim(params = set_constant_params(c(-1,beta1,gama0,0,0)), N = 100)
      data.frame(iteration = i, sim = simulated,beta1 = rep(beta1,100),gama0 = rep(gama0,100))
    })
    simulated <- do.call(bind_rows,simulated)
    heatmap.data.list[[i]] <- simulated
    i <- i + 1
    gc()
  }
}
#I wonder if this can be better written

# plan for rest of analysis:
# for each of the randomizations, I will estimate the 5 parameters using Nelder-Mead method in Optim,
# and obtain CI for each of the parameters. This will result in 100*221 estimates (whoa...). For each
# of the 221 "cells" of beta1*gama0 combinations I will have percentage of estimates for which beta1 is 
# within the CI. We will be able to examine how percentages of correct estimations chages along the
# beta1 * gama0 two-dimensional space.

# testing:
# heatmap.data.list[[1]] %>% group_by(iteration) %>% nest() %>% 
#   filter(iteration == 1) %>% 
#   mutate(estimate = lapply(data, function(x) optim(fn = count_log_like,
#                                                 par = set_params_to_optimize(numeric_vector = c(-1,mean(x$`beta1`),mean(x$`gama0`),0,0)),
#                                                 first_record_data = x$sim,const = c(0,9),hessian = T))) %>% .$estimate

clusters <- parallel::detectCores() %>%
  parallel::makeCluster()

doParallel::registerDoParallel(clusters)
parallel::clusterExport(clusters, varlist = c("count_m","count_p","count_pi","count_log_like","set_params_to_optimize"), envir = environment())
parallel::clusterCall(clusters, function() library(tidyverse))

heatmap.results.list <- lapply(heatmap.data.list,function(y){
  y %>% group_by(iteration) %>% nest() %>%  
    mutate(estimate = pbapply::pblapply(data, function(x) try(optim(fn = count_log_like,
                                                     par = set_params_to_optimize(numeric_vector = c(-1,mean(x$`beta1`),mean(x$`gama0`),0,0)),
                                                     first_record_data = x$sim,const = c(0,9),hessian = T),silent = T),cl = clusters))
})

parallel::stopCluster(clusters)


list.of.opts <- purrr::map(heatmap.results.list,3) %>% unlist(recursive = F)


list.of.params <- purrr::map(heatmap.results.list,2) %>% unlist(recursive = F)
list.of.params <- purrr::map(list.of.params, function(x) set_params_to_optimize(numeric_vector = c(-1,mean(x$`beta1`),mean(x$`gama0`),0,0)))

heatmap.success.rate <- purrr::map2(list.of.opts,list.of.params, function(x,y) try(did_optim_succeed(x,y),silent = T))

correct.rate <- map(heatmap.success.rate,6) %>% 
  map(2)
beta1.rate  <-  map(heatmap.success.rate,5) %>% 
  map(2)
gama0.rate <-  map(heatmap.success.rate,5) %>% 
  map(3)

data.for.heatmap.plot <- tibble(beta1 = unlist(beta1.rate),
       gama0 = unlist(gama0.rate),
       correct =unlist(correct.rate)) %>% 
  group_by(beta1,gama0) %>% 
  summarize(proportion = sum(correct,na.rm = T)/n())



data.for.heatmap.plot %<>% 
  mutate(level = cut(proportion,breaks = c(0,0.5,0.80,0.90,.95,1),labels = c("not significant","0.5-0.8","0.8-0.9","0.9-0.95","significant"))) %>% 
  mutate(signif = proportion > 0.95)

levels(data.for.heatmap.plot$level)

colors <- c("#177A89","#C5B4A2","#DDA35A","#A87130","#593716") 
#https://cdn1.thr.com/sites/default/files/imagecache/scale_crop_768_433/2016/02/mmfr-trl-87423.jpg
#https://www.canva.com/colors/color-palette-generator/


ggplot(data.for.heatmap.plot)+
  aes(x = beta1,y = gama0)+
  geom_tile(aes(fill = level))+
  scale_x_continuous(breaks = unique(data.for.heatmap.plot$beta1))+
  scale_y_continuous(breaks = unique(data.for.heatmap.plot$gama0))+
  scale_fill_manual(values = colors)

data.for.heatmap.plot %>% print(n=Inf)

tibble(beta1 = unlist(beta1.rate),
       gama0 = unlist(gama0.rate),
       correct =unlist(correct.rate)) %>% 
  group_by(gama0) %>% 
  summarize(proportion = sum(correct,na.rm = T)/n()) %>% {
    ggplot(.)+
      aes(x=gama0,y=proportion)+
      geom_line()
    }

did_optim_succeed_high_low <- function(opt,params){
  simulated.params <- as.data.frame(params) %>% 
    rownames_to_column()
  colnames(simulated.params) <- c("par","true")
  
  out <- t(get_parameter_se(opt)) %>% 
    as_tibble() %>% 
    mutate(par = names(opt$par)) %>% 
    left_join(simulated.params,by = "par") %>% 
    mutate(correct = (true <= `95% Upper bound` & 
                        true >= `95% Lower bound`),
           se = (`95% Upper bound` - ML)/1.96,
           estimate_bias = (true - ML)/se) %>% 
    select(par,true,estimate_bias)
  return(out)
}


heatmap.high.low <- purrr::map2(list.of.opts,list.of.params, function(x,y) try(did_optim_succeed_high_low(x,y),silent = T))

beta0.bias <- purrr::map(heatmap.high.low,3) %>% map(1)
beta1.bias <- purrr::map(heatmap.high.low,3) %>% map(2)
gama0.bias <- purrr::map(heatmap.high.low,3) %>% map(3)
gama1.bias <- purrr::map(heatmap.high.low,3) %>% map(4)
gama2.bias <- purrr::map(heatmap.high.low,3) %>% map(5)

beta1.true <- purrr::map(heatmap.high.low,2) %>% map(2)
gama0.true <- purrr::map(heatmap.high.low,2) %>% map(3)

data.for.bias.map <- tibble(beta1 = unlist(beta1.true),
                            gama0 = unlist(gama0.true),
                            beta0.bias = unlist(beta0.bias),
                            beta1.bias = unlist(beta1.bias),
                            gama0.bias = unlist(gama0.bias),
                            gama1.bias = unlist(gama1.bias),
                            gama2.bias = unlist(gama2.bias))

data.for.bias.map %>% 
  #mutate(level = cut(beta1.bias,breaks = c(-Inf,-5,-1.96,1.96,5,Inf),
   #                  labels = c("Very Under","Under","Correct","Over","Very Over"))) %>% 
  group_by(beta1,gama0) %>% nest() %>% 
  mutate(correlation = map(data,function(x) cor.test(x$gama2.bias,x$beta1.bias)[[4]])) %>% 
  unnest(correlation) %>% unnest(data) %>% 
         ggplot()+
  aes(x = beta1,y = gama0)+
  geom_tile(aes(fill = correlation))


data.for.bias.map %>% 
  group_by(beta1,gama0) %>% 
  summarize(mean.bias = mean(beta1.bias,na.rm = T)) %>% 
  mutate(bias = cut(mean.bias,breaks = c(-Inf,-1.96,1.96,Inf),
                    labels = c("Under","Correct","Over"))) %>% 
  ggplot()+
  aes(x = beta1,y = gama0)+
  geom_tile(aes(fill = bias))


data.for.bias.map %>% 
  group_by(beta1,gama0) %>% 
  summarize(under = sum(gama0.bias < -1.96,na.rm = T)/100) %>% 
  ggplot()+
  aes(x=beta1,y = gama0)+
  geom_tile(aes(fill = under))


data.for.bias.map %>% 
  group_by(beta1,gama0) %>% 
  summarize(over = sum(gama0.bias > 1.96,na.rm = T)/100) %>% 
  ggplot()+
  aes(x=beta1,y = gama0)+
  geom_tile(aes(fill = over))


data.for.bias.map %>% 
  group_by(beta1,gama0) %>% 
  summarize(over = sum(gama0.bias > 1.96,na.rm = T)/100) %>% 
  ggplot()+
  aes(x=beta1,y = over)+
  geom_point()+geom_smooth(method = "lm")

data.for.bias.map %>% 
  mutate(beta = exp(beta0.bias + beta1.bias),
         gamathing = (exp(gama0.bias + gama1.bias + gama2.bias)),
         gama = gamathing/(1+gamathing)) %>% 
  group_by(beta1,gama0) %>% 
  summarize(cor = cor.test(beta,gama)[[4]]) %>% 
  ggplot()+
  aes(x=beta1,y = gama0)+
  geom_tile(aes(fill = cor))
