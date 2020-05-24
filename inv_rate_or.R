
#install.packages("pracma")
#install.packages("neldermead")
#install.packages("matrixStats")

library(pracma)
library(neldermead)
library(matrixStats)


count_pi  <- function(S,i,params){
  # This function calculates the variable pi from Solow and Costello (2004)
  
  #pi0 = params[3]
  #pi1 = params[4]
  #pi2 = params[5]
  pi0 = params[2]
  pi1 = 0
  pi2 = 0
  
  
  num = exp(pi0 + pi1*i + pi2*exp(i-S))
  f = is.infinite(num) # if the number is infinity
  pi_total = num/(1+num)
  
  pi_total[f] <- 1
  
  return(pi_total)
}


count_p <- function(i,params){
  # This function calculates the value p from Solow and Costello (2004)
  # It uses matrix coding for efficiency
  
  S <- repmat(c(1:i),i,1)
  
  thing <- 1-count_pi(S,Conj(t(S)),params)
  thing[i,] <- 1 
  up = ones(size(thing)[1],size(thing)[2]) 
  upperones = triu(up,1) # upper triangular
  thing2 = tril(thing) + upperones
  product = colProds(as.matrix(thing2))
                     
  pst = product*count_pi(c(1:i),i,params)
                  
  return(pst)
  
}


count_m <- function(i,params){
  # This function calculates the mean, mu, from Solow and Costello (2004)
  
  m0 = params[1]
  #m1 = params[2]
  m1 = 0
  m = exp(m0 + m1*i)
  
  return(m)
}


count_log_like<- function(params){
  
  
  # This function calculates the log likelihood function for Solow and Costello (2004).
  # It takes into account any possible restrictions (See below)
  
  # params is a vector of parameters
  # restrict is a vector (same size as params) that places restrictions on the parameters.  
  # If restrict(i)=99, then there is no restriction for the ith parameter.  
  # If restrict(i)=0 (for example) then the restriction is exactly that.
  
  
  
  
  #f <- which(restrict!=99)
  #g <- which(restrict==99)
  
  #params[g] = params[g]
  #params[f] = restrict[f]
  
  lambda<-c()
  summand2<-c()
  for (i in 1:length(num_discov)){
    S <- 1:i;
    Am <- count_m(S,params) 
    Ap <- count_p(i,params)
    lambda[i] <- sum(Am * Conj(t(as.matrix(Ap))))
    
    summand2[i] <- num_discov[i]*log(lambda[i]) - lambda[i]
  }
  
  LL <- (-sum(summand2)) 
  
  return(LL) # lambda 
}


count_lambda <- function(params,N){
  # This function calculates lambda from Solow and Costello, 2004.
  # params is a vector of parameters
  lambda<-vector(mode = "numeric",length = N)
  for(i in 1:N){
    S=c(1:i)
    Am = count_m(S,params) 
    Ap = count_p(i,params)
    lambda[i] = sum(Am * Conj(t(Ap)))
  }

  return(lambda)
}



# Main Code
#########################################################################################

#This code performs the main calculations for Solow and Costello (2004) 



#NumDis <- read.csv(NumDis.txt)   # load the text file of the number of discoveries by year

NumDis <- read.table("NumDis.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

NumDis <- as.vector(NumDis)

year <- c(1851:1995)   # the time period over which discoveries were made

num_discov <- Conj(t(NumDis)) # will need to be called by some function files

# print('SF data set')

#guess <- c(0, 0, -1, 0, 0) 

guess <- c(-1.1106, 0.0135, -1.4534, 0, 0)# an initial guess 

#constr <- 99*ones(size(guess)[1],size(guess)[2]) # to set a constraint make some of these different than 99

options <- optimset(TolFun = 0.01,TolX = 0.01) # set the tolerances 

# vec1_val1 <- fminbnd(fun = count_log_like, 
#                      x0 = guess,
#                      xmin = c(-Inf, -Inf, -Inf, -Inf, -Inf),
#                      xmax = c(+Inf, +Inf, +Inf, +Inf, +Inf), options = options)# R routine for conducting MLE


vec1_val1 <- pracma::fminsearch(fn = count_log_like, x0 = guess,method="Nelder-Mead")
vec1_val2 <- fminsearch(fun = count_log_like, x0 = guess,options = options)
vec1_val3 <- optim(fn = count_log_like,par = guess,
                   lower=c(-Inf, -0.001, -Inf, -0.001, -0.001),
                   upper=c(Inf, 0.001, Inf, 0.001, 0.001 ),
                   method =  "L-BFGS-B")
vec1_val3 <- optim(fn = count_log_like,par = c(-1.1,-1.4))


vec1_val1$optbase$xopt
vec1_val1$xmin


############################################################################################################


#C1 <- count_lambda(vec1_val1[[1]],length(num_discov)); # Calculates the mean of Y


C1 <- count_lambda(vec1_val1$xmin,length(num_discov))
C2 <- count_lambda(vec1_val2$optbase$xopt,length(num_discov))
C3 <- count_lambda(vec1_val3$par,length(num_discov))
# Create the plot 
plot(year,cumsum(num_discov),'k-',year,cumsum(C1),'k--')

plot <- data.frame(t(rbind(year,num_discov,C3)))
library(tidyverse)

ggplot(plot)+
  geom_line(aes(x=year,y=cumsum(V1)))+
  geom_line(aes(x=year,y=cumsum(C3)))+
  theme_bw()

legend('Discoveries','Unrestricted')
xlabel('Year')
ylabel('Cumulative Discovery')


#################################################################################################################




temp <- fminbnd






