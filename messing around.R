source("required packages.R")
source("solow_costello_functions.R")
source("constant_p1_functions.R")

constant_p1_params <- set_constant_params(c(-1.5,0.05,-3,-0.001,0.0002),
                    parameters = c("beta0","beta1","gama0","gama1","gama2"))

params_p1 <- set_params_to_optimize(c(-1,0,-1,0.002,0.001),
                                    c("beta0","beta1","gama0","gama1","gama2"))

sampling_effort <- rnorm(100,1,1)

simulation <- sim(params = constant_p1_params, N = 100)
simulation <- sim_verbose(params = constant_p1_params, N = 100)

simulation_random_p <- sim_random_p(params = constant_p1_params, N = 100)



(annual_observation_probability <- get_p_component_const(N = 100,
                                                  params = set_constant_params(c(-1.5,0.05,-1.3,0,0)),
                                                  const = 0,pi1 =  rnorm(100,mean = 0.00,sd = 0.001))
)


(annual_observation_probability <- get_p_component(N = 100,
                                                         params =constant_p1_params,
                                                         const = 0)
)

ggplot()+
  theme_classic()+
  geom_line(aes(x = 1:100, y = cumsum(simulation)))+
  geom_line(aes(x = 1:100, y = cumsum(simulation_random_p)),color = "red")
                  



(estimate <- optim(fn = count_log_like,
                  par = constant_p1_params,
                  first_record_data = unlist(purrr::map(simulation,2)),
                  hessian = T,
                  const = constant)
)

cumsum(count_lambda_const(N = length(1:100),
                          params = estimate$par,
                          const = 0,
                          pi1 = sampling_effort))


ggplot()+
  theme_classic()+
  geom_line(aes(x = 1:100, y = cumsum(simulation)))+
  geom_line(aes(x = 1:100, y = cumsum(count_lambda(N = length(1:100),
                                                  params = estimate$par,
                                                  const = 0))),color = "red")
  

count_lambda_const()