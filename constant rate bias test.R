source("required packages.R")
#Rcpp::sourceCpp("solow_costello_rcpp_functions.cpp")
source("functions.R")
#source("constant_p1_functions.R")

linear <- readRDS("linearincrease")


guess.constant.rate.all <- set_params_to_optimize(c(-1,0.015,-2,0,0))


alreadybeengama0 <- data.for.heatmap.plot.bound$gama0
alreadybeenbeta0 <-data.for.heatmap.plot.bound$beta0
comb_exists <- cbind(alreadybeenbeta0,alreadybeengama0) %>% as_tibble()

indexes_of_new_ones <- lapply(linear, function(df) !any(
  (mean(df$gama0) %in% comb_exists$alreadybeengama0) & 
         (mean(df$beta0) %in% comb_exists$alreadybeenbeta0) &
         (any(
           which(comb_exists$alreadybeengama0 == mean(df$gama0)) %in% 
             which(comb_exists$alreadybeenbeta0 == mean(df$beta0)) )))) %>% 
  unlist()

line2 = linear[indexes_of_new_ones]
line2 = line2[sample(1:length(line2),40,F)]

clusters <- parallel::detectCores() %>%
  parallel::makeCluster()

doParallel::registerDoParallel(clusters)
parallel::clusterExport(clusters, 
                        varlist = c("guess.constant.rate.all",
                                    "set_params_to_optimize"), envir = environment())
parallel::clusterCall(clusters, function() {
  library(tidyverse)
  library(sacII)
  })


results <- pbapply::pblapply(line2,function(y){
  y %>% group_by(iteration) %>% nest() %>%  
    mutate(estimate = lapply(X = data, FUN = function(x) try(optim(fn = count_log_like,
                                                         par = guess.constant.rate.all,
                                                         first_record_data = x$sim,
                                                         const = 1,
                                                         hessian = T),
                                                   silent = T)))
}
,cl = clusters)

#saveRDS(results,"linearbiasresults")

parallel::stopCluster(clusters)

list.of.opts2 <- purrr::map(results,3) %>% unlist(recursive = F)


list.of.params2 <- purrr::map(results,2) %>% unlist(recursive = F)
list.of.params2 <- purrr::map(list.of.params2, function(x) set_params_to_optimize(numeric_vector = c(mean(x$`beta0`),0,mean(x$`gama0`),0,0)))

heatmap.success.rate2 <- purrr::map2(list.of.opts2,list.of.params2, function(x,y) try(did_optim_succeed(x,y),silent = T))

correct.rate2 <- map(heatmap.success.rate2,6) %>% 
  map(2)
beta0.rate2  <-  map(heatmap.success.rate2,5) %>% 
  map(1)
gama0.rate2 <-  map(heatmap.success.rate2,5) %>% 
  map(3)

data.for.heatmap.plot2 <- tibble(beta0 = unlist(beta0.rate2),
                                gama0 = unlist(gama0.rate2),
                                correct =unlist(correct.rate2)) %>% 
  group_by(beta0,gama0) %>% 
  summarize(proportion = sum(correct,na.rm = T)/n())



data.for.heatmap.plot2 %<>% 
  mutate(level = cut(proportion,breaks = c(0,0.5,0.80,0.90,.95,1),labels = c("not significant","0.5-0.8","0.8-0.9","0.9-0.95","significant"))) %>% 
  mutate(signif = proportion > 0.95)

levels(data.for.heatmap.plot2$level)

colors <- c("#177A89","#C5B4A2","#DDA35A","#A87130","#593716") 
#https://cdn1.thr.com/sites/default/files/imagecache/scale_crop_768_433/2016/02/mmfr-trl-87423.jpg
#https://www.canva.com/colors/color-palette-generator/

data.for.heatmap.plot.bound <- bind_rows(data.for.heatmap.plot.bound,data.for.heatmap.plot2)

saveRDS(data.for.heatmap.plot.bound,"constant rate tiledata")

data.for.heatmap.plot.bound %>% 
  mutate(level = cut(proportion,breaks = c(0,0.5,0.80,0.90,.95,1),
                     labels = c("not significant","0.5-0.8","0.8-0.9","0.9-0.95","0.95+"))) %>% 
{ggplot(.)+
  aes(x = beta0,y = gama0)+
  geom_tile(aes(fill = level))+
  scale_x_continuous(breaks = unique(.$beta0))+
  scale_y_continuous(breaks = unique(.$gama0))+
  scale_fill_manual(values = colors)
}
