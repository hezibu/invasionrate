  


# 
# 
# optim(fn = count_log_like,
#       par = set_params_to_optimize(numeric_vector = c(-1.1106,-1.4534),parameters = c("beta0","gama0")),
#       const = set_constant_params(numeric_vector = c(0,0,0),parameters = c("beta1","gama1","gama2")), 
#       first_record_data = num_discov)
# 
# 
# # Main Code
#########################################################################################

#This code performs the main calculations for Solow and Costello (2004) 


source("solow_costello_functions.R")
NumDis <- read.table("NumDis.txt") # load the text file of the number of discoveries by year
num_discov <- t(NumDis)
num_discov <- as.vector(num_discov)
rm(NumDis)


year <- c(1851:1995)   # the time period over which discoveries were made


#guess <- c(0, 0, -1, 0, 0) 
#guess <- c(-1.1106, 0.0135, -1.4534, 0, 0)# an initial guess 

guess <- set_params_to_optimize(numeric_vector = c(-1.1106,-1.4534,0),parameters = c("beta0","gama0","gama1"))
constant <- set_constant_params(numeric_vector = c(0,0),parameters = c("beta1","gama2"))

guess_all <- set_params_to_optimize(numeric_vector = c(-1.1106,0.0135,-1.4534,0,0),parameters = c("beta0","beta1","gama0","gama1","gama2"))

vec1_val1 <- optim(fn = count_log_like,par = guess_all, first_record_data = num_discov,const = constant)
vec1_val2 <- optim(fn = count_log_like,par = guess, first_record_data = num_discov,const = constant)
optimx::optimx(fn = count_log_like,par = guess, first_record_data = num_discov, const = constant)

############################################################################################################


C1 <- count_lambda(params = vec1_val1$par,N = length(num_discov),const = constant)
C2 <- count_lambda(params = vec1_val2$par,N = length(num_discov),const = constant)


plot <- data.frame(t(rbind(year,num_discov,C1,C2)))
plot <- gather(plot,key = type,value = obs,2:4)
plot <- plot %>% 
  group_by(type) %>% arrange(year) %>% mutate(cs = cumsum(obs))

plot %>% ungroup() %>%
  mutate(type = fct_relevel(type,c("num_discov","C2","C1")))%>% 
ggplot()+
  aes(x = year, y = cs,group = type,linetype = type)+
  geom_line()+
  theme_classic()+theme(legend.position = "none")
  

# deprecated:
# plot(year,cumsum(num_discov))
# lines(year,cumsum(C3))
# lines(year,cumsum(C2))
# 
# legend('Discoveries','Unrestricted')
# xlabel('Year')
# ylabel('Cumulative Discovery')


#################################################################################################################