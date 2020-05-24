#Functions for the Solow and Costello (2004) models

count_pi_const  <- function(S,i,params,const,pi1){
  
  pi0 <- ifelse(test = is.na(params["gama0"]),yes = as.numeric(const["gama0"]),no = as.numeric(params["gama0"]))
  pi0 <- ifelse(test = is.na(params["gama0"]),yes = as.numeric(pi1[i]),no = as.numeric(params["gama0"]))
  pi2 <- ifelse(test = is.na(params["gama2"]),yes = as.numeric(const["gama2"]),no = as.numeric(params["gama2"]))
  
  num  <-  exp(pi0 + pi0 + pi2*exp(i-S))
  f <- is.infinite(num) # if the number is infinity
  pi_total <-  num/(1+num)
  
  pi_total[f] <- 1
  
  return(pi_total)
}

count_p_const <- function(i,params,const,pi1){
  # This function calculates the value p from Solow and Costello (2004)
  # It uses matrix coding for efficiency
  library(matrixStats)
  library(pracma)
  
  S <- repmat(c(1:i),i,1)
  
  thing <- 1-count_pi_const(S,Conj(t(S)),params,const)
  thing[i,] <- 1 
  up  <-  ones(size(thing)[1],size(thing)[2]) 
  upperones  <-  triu(up,1) # upper triangular
  thing2  <-  tril(thing) + upperones
  product  <-  colProds(as.matrix(thing2))
  
  pst  <-  product*count_pi_const(c(1:i),i,params,const,pi1)
  
  return(pst)
  
}

count_log_like_const<- function(first_record_data,params,const,pi1){
  lambda <- vector(mode = "numeric", length = length(first_record_data))
  summand2 <- vector(mode = "numeric", length = length(first_record_data))
  for (i in seq_along(first_record_data)){
    S <- 1:i;
    Am <- count_m(S,params,const) 
    Ap <- count_p_const(i,params,const,pi1)
    lambda[i] <- sum(Am * Conj(t(as.matrix(Ap))))
    
    summand2[i] <- first_record_data[i]*log(lambda[i]) - lambda[i]
  }
  
  LL <- (-sum(summand2)) 
  
  return(LL) # lambda 
}


count_lambda_const <- function(N,params,const,pi1){
  # This function calculates lambda from Solow and Costello, 2004.
  # params is a vector of parameters
  lambda<-vector(mode = "numeric",length = N)
  for(i in 1:N){
    S=c(1:i)
    Am = count_m(S,params,const) 
    Ap = count_p_const(i,params,const,pi1)
    lambda[i] = sum(Am * Conj(t(Ap)))
  }
  
  return(lambda)
}


get_p_component_const <- function(N,params,const,pi1){
  # This function calculates lambda from Solow and Costello, 2004.
  # params is a vector of parameters
  lambda<-vector(mode = "numeric",length = N)
  for(i in 1:N){
    S=c(1:i)
    Am = count_m(S,params,const) 
    lambda[i] = sum(count_p_const(i,params,const,pi1))
  }
  
  return(lambda)
}
