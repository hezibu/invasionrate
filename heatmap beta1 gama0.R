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

#not finished yet
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
           estimate_bias = ) %>% 
    select(par,`95% Lower bound`,ML,`95% Upper bound`,true,correct)
  return(out)
}



param.matrix <- bind_rows(lapply(map(heatmap.results.list,3)[[1]],function(x) try(get_parameter_se(x),silent = T))[
  sapply(lapply(map(heatmap.results.list,3)[[1]],function(x) try(get_parameter_se(x),silent = T)),class) != "try-error"
  ],.id = "column_label") %>% 
  group_by(column_label) %>%
  filter(row_number()==1) %>%
  ungroup() %>% 
  select(-1)

params.true <-  set_params_to_optimize(c(-1,as.numeric(colMeans(map(heatmap.results.list,2)[[1]][[1]])[2:3]),0,0))

SimDesign::bias(estimate = param.matrix,parameter = params.true)
