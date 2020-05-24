stwistdata <- read.csv("AlienSpecies_MultipleDBs_Masterfile_vs1.2.csv", sep=";", stringsAsFactors=FALSE)

stwist.aves <- stwistdata %>% filter(LifeForm == "Birds")

aves.first.records <- stwist.aves %>%
  group_by(SpeciesGBIF,FirstRecord) %>% 
  summarize(n=n()) %>% 
  arrange(SpeciesGBIF,FirstRecord) %>%
  top_n(wt = FirstRecord,1) %>% ungroup() %>% select(SpeciesGBIF,FirstRecord) %>% group_by(FirstRecord) %>% summarise(NumberOfNewSpecies = n())

data.frame(FirstRecord = range(aves.first.records$FirstRecord)[1]:range(aves.first.records$FirstRecord)[2]) %>% 
  left_join(aves.first.records) %>% 
  filter(FirstRecord >= 1700) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  mutate(cs = cumsum(NumberOfNewSpecies)) %>% 
  ggplot()+
  aes(x=FirstRecord,cs)+
  geom_line()

aves.timeseries <- data.frame(FirstRecord = range(aves.first.records$FirstRecord)[1]:range(aves.first.records$FirstRecord)[2]) %>% 
  left_join(aves.first.records) %>% 
  filter(FirstRecord >= 1690) %>% 
  filter(FirstRecord >= 2005) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  .$NumberOfNewSpecies


aves.params <- set_params_to_optimize(c(0,0.02,-2,0,0))


optim(fn = count_log_like,par = aves.params,
      first_record_data = aves.timeseries, const = c(0),hessian = T)

aves.model.output <- set_constant_params(c(-4.67435808,  0.02370555, -2.13642320,  5.46325997, -0.99298235))

count_lambda(N = length(1690:range(aves.first.records$FirstRecord)[2]),params = aves.model.output)

model.data <- tibble(FirstRecord =1690:range(aves.first.records$FirstRecord)[2],
                     observed = aves.timeseries,
                     model.with.gama = count_lambda(N = length(1690:range(aves.first.records$FirstRecord)[2]),params = aves.model.output),
                     model.predict = exp(aves.model.output["beta0"]+aves.model.output["beta1"]*seq_along(FirstRecord)))

get_p_component(N = length(1690:range(aves.first.records$FirstRecord)[2]),params = aves.model.output)

aves.plot.data <- model.data %>% gather(2:4,key = group,value = num) %>% 
  group_by(group) %>% arrange(FirstRecord) %>% mutate(cs = cumsum(num))

ggplot(aves.plot.data)+
  aes(x  = FirstRecord,y = cs, group = group, linetype = group,color = group)+
  geom_line(size = 1.5)

ggplot(aves.plot.data)+
  aes(x  = FirstRecord,y = num, group = group, linetype = group,color = group)+
  geom_point()

model.data %>% print(n=Inf)


##############
aves.timeseries.70s <- data.frame(FirstRecord = range(aves.first.records$FirstRecord)[1]:range(aves.first.records$FirstRecord)[2]) %>% 
  left_join(aves.first.records) %>% 
  filter(FirstRecord >= 1970) %>% 
  mutate(NumberOfNewSpecies = ifelse(is.na(NumberOfNewSpecies),0,NumberOfNewSpecies)) %>% 
  .$NumberOfNewSpecies

optims70s <- optim(fn = count_log_like,par = aves.params,
                   first_record_data = aves.timeseries.70s, const = c(0),hessian = T)

model.data.70s <- tibble(FirstRecord =1970:range(aves.first.records$FirstRecord)[2],
                     observed = aves.timeseries.70s,
                     model.with.gama = count_lambda(N = length(1970:range(aves.first.records$FirstRecord)[2]),params = optims70s$par),
                     model.predict = exp(optims70s$par["beta0"]+optims70s$par["beta1"]*seq_along(FirstRecord)))

get_p_component(N = length(1690:range(aves.first.records$FirstRecord)[2]),params = aves.model.output)

aves.plot.data.70s <- model.data.70s %>% gather(2:4,key = group,value = num) %>% 
  group_by(group) %>% arrange(FirstRecord) %>% mutate(cs = cumsum(num))

ggplot(aves.plot.data.70s)+
  aes(x  = FirstRecord,y = cs, group = group, linetype = group,color = group)+
  geom_line(size = 1.2)

ggplot(aves.plot.data.70s)+
  aes(x  = FirstRecord,y = log(num), group = group, linetype = group,color = group)+
  geom_point()

get_p_component(N = length(1970:range(aves.first.records$FirstRecord)[2]),params = optims70s$par)

