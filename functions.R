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