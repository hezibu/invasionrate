library(tidyverse)

sim.data <- sim.discovery.data %>% 
  select(n.Inv) %>% 
  mutate(yearly  = n.Inv-lag(n.Inv)) %>% 
  select(yearly) %>% 
  filter(!is.na(yearly)) %>% 
  .$yearly

source("solow_costello_functions.R")


guess_all <- set_params_to_optimize(numeric_vector = c(-1.1106,0.0135,-1.4534,0,0),parameters = c("beta0","beta1","gama0","gama1","gama2"))
guess_all <- set_params_to_optimize(numeric_vector = c(0.2,0.05,1,1,1),parameters = c("beta0","beta1","gama0","gama1","gama2"))


parameters <- optimx::optimx(fn = count_log_like,par = guess_all, first_record_data = sim.data, const = constant)