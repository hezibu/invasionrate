guess.constant.rate.all <- set_params_to_optimize(c(-1,0.015,-2,0,0))
assume.constant <- guess.constant.rate.all[-2]

(estimate <- optim(fn = count_log_like,
                   par = guess.constant.rate.all,
                   first_record_data = sim.matrix.const.rate[,2],
                   hessian = T,
                   const = constant)
)

(estimate_constant_rate <- optim(fn = count_log_like,
                   par = assume.constant,
                   first_record_data = sim.matrix.const.rate[,2],
                   hessian = T,
                   const = set_constant_params(0,"beta1"))
)

get_parameter_se(estimate)
sim.parameters.constant.rate

ggplot()+
  theme_classic()+
  geom_line(aes(x = 1:100, y=  cumsum(sim.matrix.const.rate[,2])),color = "black")+
  geom_line(aes(x = 1:100, y = cumsum(count_lambda(100,estimate$par,const))), color = "green")+
  geom_line(aes(x = 1:100, y = cumsum(count_lambda(100,estimate_constant_rate$par,
                                                   set_constant_params(0,"beta1")))),color = "red")+
  geom_line(aes(x = 1:100, y = cumsum(count_lambda(100,sim.parameters.constant.rate))),color = "blue")
