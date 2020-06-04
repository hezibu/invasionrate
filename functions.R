sim <- function(N,params,const){
  # This function calculates lambda from Solow and Costello, 2004.
  # params is a vector of parameters
  lambda<-vector(mode = "numeric",length = N)
  for(i in 1:N){
    S=c(1:i)
    Am = count_m(S,params,const) 
    Ap = count_p(i,params,const)
    Yt = rbinom(Am,1,Ap)
    lambda[i] = round(sum(Am * Conj(t(Yt))))
  }
  
  return(lambda)
}



set_params_to_optimize <- function(numeric_vector,parameters = c("beta0","beta1","gama0","gama1","gama2")){
  return(set_names(x = numeric_vector,nm = parameters))
}

set_constant_params <- function(numeric_vector,parameters = c("beta0","beta1","gama0","gama1","gama2")){
  return(set_names(x = numeric_vector,nm = parameters))
}


take_nested_timespan <- function(matrix,timespan){ # matrix of [timeseries , simulations]
  return(tail(matrix,timespan))
}


get_parameters_estimate <- function(matrix,params){ # matrix of [timeseries , simulations]
  #this function only returns the estimate and none of the performance value such as convergence code etc.
  library(purrr)
  estimates <- apply(matrix,2,function(x)
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
    gather(key = parameter,value = estimate) %>% 
    mutate(timespan = timespan) %>% 
    return()
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
