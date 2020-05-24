#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
double count_mc(int i, NumericVector params, NumericVector constant){
  double m0 = params.containsElementNamed("beta0") ? params["beta0"] : constant["beta0"]; //containsElementNamed is good for our purposes
  double m1 = params.containsElementNamed("beta1") ? params["beta1"] : constant["beta1"];
  double m  = exp(m0 + m1*i);
  return(m);
}

// [[Rcpp::export]]
double count_p(int i,NumericVector params,NumericVector constant){
 NumericMatrix a = NumericMatrix::create( )
  NumericMatix S = (c(1:i),i,1)
   
   thing = 1-count_pi(S,Conj(t(S)),params,constant)
   thing[i,] <- 1 
   up  <-  ones(size(thing)[1],size(thing)[2]) 
   upperones  <-  triu(up,1) # upper triangular
   thing2  <-  tril(thing) + upperones
   product  <-  colProds(as.matrix(thing2))
   
   pst  <-  product*count_pi(c(1:i),i,params,constant)
   
   return(pst)
   
 }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

#Functions for the Solow and Costello (2004) models

# count_pi  <- function(S,i,params,const){
#   
#   pi0 <- ifelse(test = is.na(params["gama0"]),yes = as.numeric(const["gama0"]),no = as.numeric(params["gama0"]))
#   pi1 <- ifelse(test = is.na(params["gama1"]),yes = as.numeric(const["gama1"]),no = as.numeric(params["gama1"]))
#   pi2 <- ifelse(test = is.na(params["gama2"]),yes = as.numeric(const["gama2"]),no = as.numeric(params["gama2"]))
#   
#   num  <-  exp(pi0 + pi1*i + pi2*exp(i-S))
#   f <- is.infinite(num) # if the number is infinity
#   pi_total <-  num/(1+num)
#   
#   pi_total[f] <- 1
#   
#   return(pi_total)
# }
# 
# count_p <- function(i,params,const){
# # This function calculates the value p from Solow and Costello (2004)
# # It uses matrix coding for efficiency
#   library(matrixStats)
#   library(pracma)
#   
#   S <- repmat(c(1:i),i,1)
#   
#   thing <- 1-count_pi(S,Conj(t(S)),params,const)
#   thing[i,] <- 1 
#   up  <-  ones(size(thing)[1],size(thing)[2]) 
#   upperones  <-  triu(up,1) # upper triangular
#   thing2  <-  tril(thing) + upperones
#   product  <-  colProds(as.matrix(thing2))
#   
#   pst  <-  product*count_pi(c(1:i),i,params,const)
#   
#   return(pst)
#   
# }
# 

# 
# count_log_like<- function(first_record_data,params,const){
#   lambda<-vector(mode = "numeric", length = length(first_record_data))
#   summand2<-vector(mode = "numeric", length = length(first_record_data))
#   for (i in seq_along(first_record_data)){
#     S <- 1:i;
#     Am <- count_m(S,params,const) 
#       Ap <- count_p(i,params,const)
#       lambda[i] <- sum(Am * Conj(t(as.matrix(Ap))))
#       
#       summand2[i] <- first_record_data[i]*log(lambda[i]) - lambda[i]
#   }
#   
#   LL <- (-sum(summand2)) 
#     
#     return(LL) # lambda 
# }
# 
# 
# count_lambda <- function(N,params,const){
# # This function calculates lambda from Solow and Costello, 2004.
# # params is a vector of parameters
#   lambda<-vector(mode = "numeric",length = N)
#   for(i in 1:N){
#     S=c(1:i)
#     Am = count_m(S,params,const) 
#     Ap = count_p(i,params,const)
#     lambda[i] = sum(Am * Conj(t(Ap)))
#   }
#   
#   return(lambda)
# }
# 
# sim <- function(N,params,const){
# # This function calculates lambda from Solow and Costello, 2004.
# # params is a vector of parameters
#   lambda<-vector(mode = "numeric",length = N)
#   for(i in 1:N){
#     S=c(1:i)
#     Am = count_m(S,params,const) 
#     Ap = count_p(i,params,const)
#     Yt = rbinom(Am,1,Ap)
#     lambda[i] = round(sum(Am * Conj(t(Yt))))
#   }
#   
#   return(lambda)
# }
# 
# 
# 
# set_params_to_optimize <- function(numeric_vector,parameters = c("beta0","beta1","gama0","gama1","gama2")){
#   return(set_names(x = numeric_vector,nm = parameters))
# }
# 
# set_constant_params <- function(numeric_vector,parameters = c("beta0","beta1","gama0","gama1","gama2")){
#   return(set_names(x = numeric_vector,nm = parameters))
# }
# 


*/
