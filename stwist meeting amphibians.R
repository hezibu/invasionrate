stwist.amphibians <-stwistdata %>% filter(LifeForm == "Amphibians")

amphibians.first.records <- stwist.amphibians %>%
  group_by(SpeciesGBIF,FirstRecord) %>% 
  summarize(n=n()) %>% 
  arrange(SpeciesGBIF,FirstRecord) %>%
  top_n(wt = FirstRecord,1) %>% ungroup() %>% select(SpeciesGBIF,FirstRecord) %>% group_by(FirstRecord) %>% summarise(NumberOfNewSpecies = n())

data.frame(FirstRecord = range(amphibians.first.records$FirstRecord)[1]:range(amphibians.first.records$FirstRecord)[2]) %>% 
  left_join(amphibians.first.records) %>% 
  #filter(FirstRecord >= 1700) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  mutate(cs = cumsum(NumberOfNewSpecies)) %>% 
  ggplot()+
  aes(x=FirstRecord,cs)+
  geom_line()

amphibians.timeseries <- data.frame(FirstRecord = range(amphibians.first.records$FirstRecord)[1]:range(amphibians.first.records$FirstRecord)[2]) %>% 
  left_join(amphibians.first.records) %>% 
  #filter(FirstRecord >= 1690) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  .$NumberOfNewSpecies


amphibians.params <- set_params_to_optimize(c(1,0.01,-1,0,0))


amphibians.model <- optim(fn = count_log_like,par = amphibians.params,
      first_record_data = amphibians.timeseries, const = c(0),hessian = T)


count_lambda(N = length(1690:range(amphibians.first.records$FirstRecord)[2]),params = amphibians.model.output)

model.data.amphibians <- tibble(FirstRecord =range(amphibians.first.records$FirstRecord)[1]:range(amphibians.first.records$FirstRecord)[2],
                     observed = amphibians.timeseries,
                     model.with.gama = count_lambda(N = length(range(amphibians.first.records$FirstRecord)[1]:range(amphibians.first.records$FirstRecord)[2]),params = amphibians.model$par),
                     model.predict = exp(amphibians.model$par["beta0"]+amphibians.model$par["beta1"]*seq_along(FirstRecord)))

get_p_component(N = length(1690:range(amphibians.first.records$FirstRecord)[2]),params = amphibians.model.output)

amphibians.plot.data <- model.data.amphibians %>% gather(2:4,key = group,value = num) %>% 
  group_by(group) %>% arrange(FirstRecord) %>% mutate(cs = cumsum(num))

ggplot(amphibians.plot.data)+
  aes(x  = FirstRecord,y = cs, group = group, linetype = group,color = group)+
  geom_line(size = 1.3)


ggplot(amphibians.plot.data)+
  aes(x  = FirstRecord,y = num, group = group, linetype = group,color = group)+
  geom_line(size = 1.3)


model.observed <- lm(formula = log(num+1) ~ seq_along(FirstRecord) ,data = amphibians.plot.data,subset = group == "observed")
summary(model.observed)

model.observed <- lm(formula = log(num+1) ~ seq_along(FirstRecord) ,data = amphibians.plot.data,subset = group == "observed")
summary(model.observed)


#######

amphibians.timeseries.70s <- data.frame(FirstRecord = range(amphibians.first.records$FirstRecord)[1]:range(amphibians.first.records$FirstRecord)[2]) %>% 
  left_join(amphibians.first.records) %>% 
  filter(FirstRecord >= 1970) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  .$NumberOfNewSpecies

amphibians.params <- set_params_to_optimize(c(-0.1,0.00,-1,0,0))

optims70s.amph <- optim(fn = count_log_like,par = amphibians.params,
                   first_record_data = amphibians.timeseries.70s, const = c(0),hessian = T)

model.data.70s <- tibble(FirstRecord =1970:range(amphibians.first.records$FirstRecord)[2],
                         observed = amphibians.timeseries.70s,
                         model.with.gama = count_lambda(N = length(1970:range(amphibians.first.records$FirstRecord)[2]),params = optims70s$par),
                         model.predict = exp(optims70s$par["beta0"]+optims70s$par["beta1"]*seq_along(FirstRecord)))

get_p_component(N = length(1690:range(amphibians.first.records$FirstRecord)[2]),params = amphibians.model.output)

amphibians.plot.data.70s <- model.data.70s %>% gather(2:4,key = group,value = num) %>% 
  group_by(group) %>% arrange(FirstRecord) %>% mutate(cs = cumsum(num))

ggplot(amphibians.plot.data.70s)+
  aes(x  = FirstRecord,y = cs, group = group, linetype = group,color = group)+
  geom_line()

