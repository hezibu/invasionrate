source("required packages.R")
source("solow_costello_functions.R")

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
        first_record_data = taxon.timeseries, const = c(0),hessian = T,method = "BFGS")
  
  param.samps <- MASS::mvrnorm(n = 1000,mu = taxon.model.output$par,Sigma = solve(taxon.model.output$hessian),tol = 0.0001)
  predictions <- apply(param.samps,1,function(par) cumsum(count_lambda(N = length(timeframe[1]:timeframe[2]),params = par,const = 0)))


  predictions.sd <- apply(predictions,1,sd)
  
  
  exponent.predictions <- apply(param.samps[,1:2],1,function(row) 
    exp(row["beta0"]+row["beta1"]*(seq_along(timeframe[1]:timeframe[2]))))
  
  exp.predictions.sd <- apply(exponent.predictions,1,sd)
  
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
    group_by(group) %>% arrange(first.record) %>% mutate(cs = cumsum(num)) %>% dplyr::ungroup() 
  
  plot1 <- taxon.plot.data %>% 
    select(first.record,group,num) %>% 
    filter(group %in% c("simple.exponent",'snc.exponent')) %>%
    ggplot()+
    geom_line(aes(x  = first.record,y = num,linetype =group), size = 1.5)+
    geom_ribbon(data = subset(test,group == "snc.exponent"),aes(x = first.record, 
                                                                ymax = num + 1.96 * exp.predictions.sd,
                                                                ymin = num - 1.96 * exp.predictions.sd),alpha = 0.2)+
    xlab("Year")+ylab("Cummulative Number of Species")+
    scale_linetype_discrete(labels = c("Observed","Model Fit"))+
    theme_classic()+
    theme(legend.title = element_blank(),legend.position = c(.15,.85))
  
  plot2 <- taxon.plot.data %>% 
    select(first.record,group,cs) %>% 
    filter(group %in% c("observed",'snc.fit')) %>% 
    ggplot()+
    geom_line(aes(x  = first.record,y = cs,linetype =group), size = 1.5)+
    geom_ribbon(data = subset(test,group == "snc.fit"),aes(x = first.record,
                                                           ymax = cs + 1.96 * predictions.sd,
                                                           ymin = cs - 1.96 * predictions.sd),alpha = 0.2)+
    xlab("Year")+ylab("Cummulative Number of Species")+
    scale_linetype_discrete(labels = c("Observed","Model Fit"))+
    theme_classic()+
    theme(legend.title = element_blank(),legend.position = c(.15,.85))
    
  
  return(list(taxon.model.output,plot1,plot2))
}


birds <- get_plot_stwist("Birds",taxon.params = set_params_to_optimize(c(1,0.01,-1,0,0)),
                timeframe = c(1800,2000))

