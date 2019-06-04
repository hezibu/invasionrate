first.record.hanno <- read.csv("hannodata.csv")

first.record.hanno %>% 
 
  group_by(FirstRecord) %>% 
  summarize(number_of_new_species = n()) %>% 
  {ggplot(.)+
      aes(x=FirstRecord,y = cumsum(number_of_new_species))+
      geom_point()}


new.records.per.region <- first.record.hanno %>% 
  filter(FirstRecord > 1500) %>% 
  group_by(Region,FirstRecord) %>% 
  summarize(number_of_new_species = n())

hanno.plot <- new.records.per.region %>% group_by(Region) %>% nest() %>% 
  mutate(plot = map2(data,Region, function(x,y) {
    ggplot(x) + aes(x=FirstRecord, y=cumsum(number_of_new_species))+
      geom_point()+ggtitle({y})
  }))

new.records.plot <- gridExtra::arrangeGrob(grobs = hanno.plot$plot,ncol = 10)
ggsave(new.records.plot,filename = "new_records.png",
       width = 50, height =  80,limitsize = F)


australia <- new.records.per.region %>% filter(Region == "Australia")

australia.first.records <- data.frame(FirstRecord = range(australia$FirstRecord)[1]:range(australia$FirstRecord)[2]) %>% 
  left_join(australia,by ="FirstRecord") %>% 
  select(-Region) %>% 
  mutate(number_of_new_species = ifelse(is.na(number_of_new_species),0,number_of_new_species)) %>% 
  .$number_of_new_species


plot(1:length(australia.first.records),cumsum(australia.first.records))


australia.params <- set_params_to_optimize(c(-1,0.006,-1.4,0,0))

austalia.estimate <- optim(count_log_like,australia.params,first_record_data = australia.first.records,constant = australia.params,hessian = T)


get_parameters_estimate(t(repmat(australia.first.records,2,1)),params = australia.params,const = australia.params)

australia.first.records <- as.integer(australia.first.records)

optimx::optimx(par = australia.params,fn = count_log_like,hess = T,first_record_data = australia.first.records,const = australia.params)

# c(1.008250, 0.010535342,  1.298154, 0.9463893,  2.327418)
# 
# set_params_to_optimize(c(1.008250, 0.010535342,  1.298154, 0.9463893,  2.327418)
# )

australia.predict <- count_lambda(params = set_params_to_optimize(c(1.008250, 0.010535342,  1.298154, 0.9463893,  2.327418)
),N = length(australia.first.records),const = constant)

aus.plot <- data.frame(t(rbind( range(australia$FirstRecord)[1]:range(australia$FirstRecord)[2],
                                australia.first.records,
                                australia.predict)))
aus.plot <- gather(aus.plot,key = type,value = obs,2:3)
aus.plot <- aus.plot %>% 
  group_by(type) %>% arrange(V1) %>% mutate(cs = cumsum(obs))

aus.plot %>% ungroup() %>%
  #mutate(type = fct_relevel(type,c("num_discov","C2","C1")))%>% 
  ggplot()+
  aes(x = V1, y = cs,group = type,linetype = type)+
  geom_line()+
  theme_classic()+theme(legend.position = "none")


#
optimx::optimx(par = australia.params[c(1,3)],fn = count_log_like,first_record_data = australia.first.records,const = set_constant_params(c(0,0,0,0,0)))

australia.params.2 <- set_params_to_optimize(c(-1,-0.05,-2.5,0,0))

optimx::optimx(par = australia.params.2,fn = count_log_like,first_record_data = australia.first.records,const = set_constant_params(c(0,0,0,0,0)))
c(1.054028, 0.01072071, -1.9761064, -1.400645, 0.6854810)
australia.predict.2 <- count_lambda(params = set_params_to_optimize(c(1.054028, 0.01072071, -1.9761064, -1.400645, 0.6854810)
),N = length(australia.first.records),const = constant)


aus.plot <- data.frame(t(rbind( year = range(australia$FirstRecord)[1]:range(australia$FirstRecord)[2],
                                australia.first.records,
                                australia.predict,australia.predict.2)))
aus.plot <- gather(aus.plot,key = type,value = obs,2:4)
aus.plot <- aus.plot %>% 
  group_by(type) %>% arrange(year) %>% mutate(cs = cumsum(obs))

aus.plot %>% ungroup() %>%
  #mutate(type = fct_relevel(type,c("num_discov","C2","C1")))%>% 
  ggplot()+
  aes(x = year, y = cs,group = type,linetype = type)+
  geom_line()+
  theme_classic()+theme(legend.position = "none")

# c(1.008250, 0.010535342,  1.298154, 0.9463893,  2.327418)
# c(1.054028, 0.01072071, -1.9761064, -1.400645, 0.6854810)
