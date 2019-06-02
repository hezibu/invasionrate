sim.parameters.timespan = set_constant_params(c(-1.1,0.014,-1.46,0.00001,0.0000004),
                                     parameters = c("beta0","beta1","gama0","gama1","gama2"))

sim.matrix <- replicate(expr = sim(params = sim.parameters.timespan, N = 150),n = 100)

library(pacman)
p_load(pbapply,optimx,purrr)

paramateres.estimate.10.years <- pbapply(sim.matrix.10.years,2,function(x)
        as.numeric(optim(fn = count_log_like,
                         par = sim.parameters.timespan,
                         first_record_data = x,
                         const = constant)$par))


take_nested_timespan <- function(matrix,timespan){ # matrix of [timeseries , simulations]
  return(tail(matrix,timespan))
}


get_parameters_estimate_only_par <- function(matrix,params){ # matrix of [timeseries , simulations]
  #this function only returns the estimate and none of the performance value such as convergence code etc.
  library(purrr)
  estimates <- pbapply(matrix,2,function(x)
    as.numeric(optim(fn = count_log_like,
                     par = params,
                     first_record_data = x,
                     const = constant)$par))
  estimates <- t(estimates) %>% 
    as_data_frame() %>%
    set_names(nm = names(params))
  return(estimates)
}

get_parameters_specified_timespan <- function(matrix,params,timespan){
  take_nested_timespan(matrix,timespan) %>% 
    get_parameters_estimate(params) %>% 
    return()
}

get_parameter_se <- function(opt){
  #this method is based on Nathaniel Philips article located here:
  #https://rstudio-pubs-static.s3.amazonaws.com/107801_2785d7d7a49744539ef21eaaebe4fe5a.html
  hessian <- opt$hessian 
  hessian.inv <- solve(hessian)
  p.se <- sqrt(diag(hessian.inv))
  CI.matrix <- as.data.frame(matrix(NA, nrow = 3, ncol = length(opt$par)))
  
  CI.matrix[1,] <- opt$par
  CI.matrix[2,] <- opt$par - 1.96 * p.se
  CI.matrix[3,] <- opt$par + 1.96 * p.se
  names(CI.matrix) <- names(opt$par)
  rownames(CI.matrix) <- c("ML", "95% Lower bound", "95% Upper bound")
  
  CI.matrix
}

parameters.results.timespans <- list()
timespan.intervals <- seq(10,150,by=10)
for (i in seq_along(timespan.intervals)){
  print(i)
  parameters.results.timespans[[i]] <- get_parameters_specified_timespan(sim.matrix,
                                                                         sim.parameters.timespan,
                                                                         timespan.intervals[i])
}


clusters <- 7 %>%
  parallel::makeCluster()
doParallel::registerDoParallel(clusters)
parallel::clusterExport(clusters, varlist = c("set_params_to_optimize","sim.matrix","sim.parameters.timespan","get_parameters_estimate","take_nested_timespan","get_parameters_specified_timespan","count_m","count_p","count_pi","count_log_like","set_params_to_optimize"), envir = environment())
parallel::clusterCall(clusters, function() library(tidyverse))
parameters.results.timespans <- pblapply(timespan.intervals,
  function(i) get_parameters_specified_timespan(sim.matrix,set_params_to_optimize(c(-1,0,-1.5,0,0)), i),cl = clusters)
parallel::stopCluster(clusters)

# true.values <- as.data.frame(sim.parameters.timespan) %>%
#   mutate(parameter = row.names(.)) %>%
#   rename(true_estimate = sim.parameters.timespan) %>% 
#   select(parameter,true_estimate)


timespan.success.rate <- lapply(parameters.results.timespans, function(x) {
  lapply(x,function(y) try(did_optim_succeed(y,sim.parameters.timespan),silent = T))
})



timespan.success.rate.merged.dataframes <- 
  lapply(X = timespan.success.rate, FUN = function(x) bind_rows(x[sapply(x,class) != "try-error"]))


timespan.success.rate <- lapply(timespan.success.rate.merged.dataframes, function(x){
  x %>% 
    group_by(par) %>% 
    summarise(prop = sum(correct,na.rm = T)/100)
})

timespan.success.rate <- lapply(seq_along(timespan.intervals), function(i) {timespan.success.rate[[i]] %>% 
    mutate(timespan = timespan.intervals[[i]])})

timespan.success.rate <- do.call(bind_rows,timespan.success.rate)

ggplot(timespan.success.rate) +
  aes(x=timespan,y = prop,color = par,group = par)+
  geom_line(size = 1.3) + theme_bw()+geom_hline(yintercept = 0.95,color = "red",linetype = 2,size = 1.2)


# complete.data.ts %>% group_by(timespan,parameter) %>% 
#   summarise(mean = mean(estimate),
#             c95 = quantile(estimate,0.95),
#             c05 = quantile(estimate,0.05)) %>% 
#   left_join(true.values,by = "parameter") %>% 
#   #mutate(ses = (true_estimate - mean)/sd) %>% 
#   ggplot()+
#   aes(color = parameter)+
#   geom_point(aes(x=timespan,y=mean))+theme_classic()+geom_ribbon(aes(x=timespan,ymin = c05, ymax = c95,fill = parameter,alpha = 0.2))+
#   geom_line(aes(x=timespan,y= true_estimate),size = 1.4,position = position_dodge(width = 0.2))+
#   coord_cartesian(ylim = c(-0.2,0.2))+
#   facet_grid(~parameter)
# 


################

timespan.heatmap.data.list <-  vector(mode = "list",length = length(seq(-0.01,0.02,by = 0.0010)))

i <- 1
for (beta1 in seq(-0.01,0.02,by = 0.0010)){
  timespan.simulated <- lapply(1:100, function(i) {
    timespan.simulated <- sim(params = set_constant_params(c(-1,beta1,-1,0,0)), N = 100)
    data.frame(iteration = i, sim = timespan.simulated,beta1 = rep(beta1,100))
  })
  timespan.simulated <- do.call(bind_rows,timespan.simulated)
  timespan.heatmap.data.list[[i]] <- timespan.simulated
  i <- i + 1
  gc()
}

timespan.heatmap.data.list <- 
  pbapply::pblapply(seq(-0.01,0.02,by = 0.0010),
                                                function(beta1) replicate(sim(params = set_constant_params(c(-1,beta1,-1,0,0)),
                                                                              N = 100),n = 100))

#I will later name the elements according to the beta1 value

clusters <- parallel::detectCores() %>%
  parallel::makeCluster()
doParallel::registerDoParallel(clusters)
parallel::clusterExport(clusters, varlist = c("timespan.heatmap.data.list","timespan.intervals","set_params_to_optimize","sim.matrix","sim.parameters.timespan",
                                              "get_parameters_estimate","take_nested_timespan",
                                              "get_parameters_specified_timespan","count_m","count_p","count_pi",
                                              "count_log_like","set_params_to_optimize"), envir = environment())
parallel::clusterCall(clusters, function() library(tidyverse))

timespan.heatmap.results.list <- lapply(seq_along(seq(-0.01,0.02,by = 0.0010)),
                             function(beta1) pblapply(timespan.intervals[c(1,4,6,10)],
                             function(i) try(get_parameters_specified_timespan(matrix = timespan.heatmap.data.list[[beta1]],
                             params = set_params_to_optimize(c(-1,seq(-0.01,0.02,by = 0.0010)[beta1],-1,0,0)),
                             timespan =  i),silent = T),cl = clusters))

parallel::stopCluster(clusters)

purrr::map(timespan.heatmap.results.list,1) %>% unlist(recursive = F)

timespan.heatmap.results.list %>% unlist(recursive = F)

bind_rows(timespan.heatmap.results.list)

optim(fn = count_log_like,par = set_params_to_optimize(c(-1,0.02,-1,0,0)),
      first_record_data =timespan.heatmap.data.list[[31]][,100],const = c(0))


timespan.heatmap.results.list[[1]][[1]][[2]]

heatmap.result.dataframe <- data_frame(beta1 = rep(seq(-0.01,0.02,by = 0.0010),10*100),
       timespan = rep(timespan.intervals[1:10],31*100),
       index = rep(1:100,each = 31*10)) %>% 
  arrange(beta1,timespan,index)

index <- 1
list.r <- list()
for (i in 1:31){
  for (j in 1:10){
    for (k in 1:100){
      list.r[[index]] <- try(silent = T,did_optim_succeed(params = set_params_to_optimize(c(-1,seq(-0.01,0.02,by = 0.0010)[i],-1,0,0)),
                            timespan.heatmap.results.list %>% .[[i]] %>% .[[j]] %>% .[[k]]))
      index <- index +1
    }
  }
}



list.r[[1003]]

heatmap.result.dataframe$success <- list.r

timespan.correct.rate <- map(list.r,6) %>% 
  map(2)
timespan.beta1.rate  <-  map(list.r,5) %>% 
  map(2)

heatmap.result.dataframe$beta1correct <- map(list.r,6) %>% map(2)


heatmap.result.dataframe$beta1correct[sapply(heatmap.result.dataframe$beta1correct,is.null)] <- NA


heatmap.result.dataframe <- unnest(heatmap.result.dataframe,beta1correct)



heatmap.result.dataframe %>% group_by(beta1,timespan) %>% summarize(correctrate = sum(beta1correct,na.rm = T)/100) %>% 
  print(n=Inf)

unlist(heatmap.result.dataframe$beta1correct)



