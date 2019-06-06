get_plot_stwist <- function(taxon,taxon.params,timeframe = c(1900,2000)) {
  stwistdata <- read.csv("AlienSpecies_MultipleDBs_Masterfile_vs1.2.csv", sep=";", stringsAsFactors=FALSE)
  
  taxon.stwist <- stwistdata %>% filter(LifeForm == taxon)
  
  taxon.first.records <- taxon.stwist %>%
    rename(first.record = FirstRecord) %>% 
    group_by(SpeciesGBIF,first.record) %>% 
    summarize(n=n()) %>% 
    arrange(SpeciesGBIF,first.record) %>%
    top_n(wt = first.record,1) %>% ungroup() %>% select(SpeciesGBIF,first.record) %>%
    group_by(first.record) %>% summarise(number.of.new.species = n())
  
    taxon.first.records <- taxon.first.records %>% 
      filter(first.record >= timeframe[1]) %>% 
      filter(first.record <= timeframe[2])
    
  
  taxon.timeseries.data <- data.frame(first.record = timeframe[1]:timeframe[2]) %>% 
    left_join(taxon.first.records,by = "first.record") %>% 
    mutate(number.of.new.species = ifelse(is.na(number.of.new.species),0,number.of.new.species),
           years.from.start = seq_along(first.record)-1)
  
  
  
  taxon.timeseries <- taxon.timeseries.data$number.of.new.species

  taxon.model.output <- optim(fn = count_log_like,par = taxon.params,
        first_record_data = taxon.timeseries, const = c(0),hessian = T)
  
  linear.model <- lm(log(number.of.new.species+1) ~ years.from.start,
                     data = taxon.timeseries.data)
  
  start <- list(beta0 = coef(linear.model)[1],
                beta1 = coef(linear.model)[2])

  simple.model <- nls(number.of.new.species ~ exp(beta0 + beta1*years.from.start),
                      taxon.timeseries.data,start = start)
  
  s.beta0 <- coef(simple.model)[1]
  s.beta1 <- coef(simple.model)[2]
  
  
  simple.model.predict <- exp(s.beta0 + s.beta1*taxon.timeseries.data$years.from.start)

  
  model.data <- tibble(first.record = timeframe[1]:timeframe[2],
                       observed = taxon.timeseries,
                       simple.exponent = simple.model.predict,
                       snc.fit = count_lambda(N = length(taxon.timeseries),params = taxon.model.output$par),
                       snc.exponent = exp(taxon.model.output$par["beta0"]+taxon.model.output$par["beta1"]*(seq_along(first.record)-1)))

  taxon.plot.data <- model.data %>% gather(2:5,key = group,value = num) %>%
    group_by(group) %>% arrange(first.record) %>% mutate(cs = cumsum(num)) %>% dplyr::ungroup() %>% 
    mutate(group,forcats::fct_relevel(group,c("observed","simple.exponent","snc.fit","snc.exponent")))


  plot1 <- taxon.plot.data %>% 
    filter(group %in% c("simple.exponent",'snc.exponent')) %>% 
    ggplot()+
    aes(x  = first.record,y = num, group = group, linetype = group)+
    geom_line(size = 1.5)+xlab("Year")+ylab("Yealy Number of Species")+
    scale_linetype_discrete(labels = c("Observed Rate","Model Derived Rate"))
  
  plot2 <- taxon.plot.data %>% 
    filter(group %in% c("observed",'snc.fit')) %>% 
    ggplot()+
    aes(x  = first.record,y = cs, group = group, linetype = group)+
    geom_line(size = 1.5)+xlab("Year")+ylab("Cummulative Number of Species")+
    scale_linetype_discrete(labels = c("Observed","Model Fit"))
  
  return(list(taxon.model.output,plot1,plot2))
}



get_plot_stwist("Birds",taxon.params = aves.params,
                timeframe = c(1970,2000))

