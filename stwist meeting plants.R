stwist.plant <-stwistdata %>% filter(LifeForm == "Vascular plants")

plant.first.records <- stwist.plant %>%
  group_by(SpeciesGBIF,FirstRecord) %>% 
  summarize(n=n()) %>% 
  arrange(SpeciesGBIF,FirstRecord) %>%
  top_n(wt = FirstRecord,1) %>% ungroup() %>% select(SpeciesGBIF,FirstRecord) %>% group_by(FirstRecord) %>% summarise(NumberOfNewSpecies = n())

data.frame(FirstRecord = range(plant.first.records$FirstRecord)[1]:range(plant.first.records$FirstRecord)[2]) %>% 
  left_join(plant.first.records) %>% 
  filter(FirstRecord >= 1700) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  mutate(cs = cumsum(NumberOfNewSpecies)) %>% 
  ggplot()+
  aes(x=FirstRecord,cs)+
  geom_line()

plant.timeseries <- data.frame(FirstRecord = range(plant.first.records$FirstRecord)[1]:range(plant.first.records$FirstRecord)[2]) %>% 
  left_join(plant.first.records) %>% 
  filter(FirstRecord >= 1690) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  .$NumberOfNewSpecies


plant.params <- set_params_to_optimize(c(1,0.03,-1,0,0))


plant.model <- optim(fn = count_log_like,par = plant.params,
                          first_record_data = plant.timeseries, const = c(0),hessian = T)


count_lambda(N = length(1690:range(plant.first.records$FirstRecord)[2]),params = plant.model.output)

model.data.plant <- tibble(FirstRecord =1690:range(plant.first.records$FirstRecord)[2],
                                observed = plant.timeseries,
                                model.with.gama = count_lambda(N = length(1690:range(plant.first.records$FirstRecord)[2]),params = plant.model$par),
                                model.predict = exp(plant.model$par["beta0"]+plant.model$par["beta1"]*seq_along(FirstRecord)))

get_p_component(N = length(1690:range(plant.first.records$FirstRecord)[2]),params = plant.model.output)

plant.plot.data <- model.data.plant %>% gather(2:4,key = group,value = num) %>% 
  group_by(group) %>% arrange(FirstRecord) %>% mutate(cs = cumsum(num))

ggplot(plant.plot.data)+
  aes(x  = FirstRecord,y = num, group = group, linetype = group,color = group)+
  geom_line()

model.observed <- lm(formula = log(num+1) ~ seq_along(FirstRecord) ,data = plant.plot.data,subset = group == "observed")
summary(model.observed)


############

plant.timeseries.70s <- data.frame(FirstRecord = range(plant.first.records$FirstRecord)[1]:range(plant.first.records$FirstRecord)[2]) %>% 
  left_join(plant.first.records) %>% 
  filter(FirstRecord >= 1970) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  .$NumberOfNewSpecies



optims70s <- optim(fn = count_log_like,par = plant.params,
                   first_record_data = plant.timeseries.70s, const = c(0),hessian = T)

model.data.70s <- tibble(FirstRecord =1970:range(plant.first.records$FirstRecord)[2],
                         observed = plant.timeseries.70s,
                         model.with.gama = count_lambda(N = length(1970:range(plant.first.records$FirstRecord)[2]),params = optims70s$par),
                         model.predict = exp(optims70s$par["beta0"]+optims70s$par["beta1"]*seq_along(FirstRecord)))

get_p_component(N = length(1690:range(plants.first.records$FirstRecord)[2]),params = plants.model.output)

plant.plot.data.70s <- model.data.70s %>% gather(2:4,key = group,value = num) %>% 
  group_by(group) %>% arrange(FirstRecord) %>% mutate(cs = cumsum(num))

ggplot(plant.plot.data.70s)+
  aes(x  = FirstRecord,y = num, group = group, linetype = group,color = group)+
  geom_line()


#2005 is last timestep
