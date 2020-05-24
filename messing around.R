source("required packages.R")
source("solow_costello_functions.R")
source("constant_p1_functions.R")

constant_p1_params <- set_constant_params(c(-1.1,0,0,0,0.0000004),
                    parameters = c("beta0","beta1","gama0","gama1","gama2"))

params_p1 <- set_params_to_optimize(c(-1,0.01,-1,0.002,0.001),
                                    c("beta0","beta1","gama0","gama1","gama2"))

sampling_effort <- rnorm(100,1,1)

simulation <- sim(params = constant_p1_params, N = 100)



annual_observation_probability <- get_p_component(N = 100,
                                                  params = set_constant_params(c(-1.1,0.014,0,0,0.0000004)),                                                 
                                                  const = 0)

ggplot()+
  theme_classic()+
  geom_line(aes(x = 1:100, y = cumsum(simulation)))

(estimate <- optim(fn = count_log_like,
                  par = params_p1,
                  first_record_data = simulation,
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