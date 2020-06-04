#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat triu(arma::mat m,int n){
  arma::mat Z = m;
  for(int i = 0; i < Z.n_rows; ++i){
    for(int j = 0; j < Z.n_cols; ++j){
      if((i + n - 1) >= j) Z(i,j) = 0;
    }
  }
  return(Z);
}

// [[Rcpp::export]]
arma::mat tril(arma::mat m,int n){
  arma::mat Z = m;
  for(int i = 0; i < Z.n_rows; ++i){
    for(int j = 0; j < Z.n_cols; ++j){
      if((i + n) < j) Z(i,j) = 0;
    }
  }
  return(Z);
}

// [[Rcpp::export]]
arma::rowvec colProds(arma::mat m){
  arma::rowvec a = ones(m.n_cols).t();
  for(int i = 0; i < m.n_cols; ++i){
    for(int j = 0; j < m.n_rows; ++j){
      a.col(j) = a.col(j) *  m(i,j);
    }
  }
  return a;
}

// [[Rcpp::export]]
arma::mat count_pi_mat(arma::mat S, arma::mat i, NumericVector params, NumericVector constant) {
  
  double pi0 = params.containsElementNamed("gama0") ? params["gama0"] : constant["gama0"]; //containsElementNamed is good for our purposes
  double pi1 = params.containsElementNamed("gama1") ? params["gama1"] : constant["gama1"];
  double pi2 = params.containsElementNamed("gama2") ? params["gama2"] : constant["gama2"];
  
  
  arma::mat num =  exp(pi0 + pi1*i + pi2*exp(i - S));
  arma::mat pi_total =  num/(1+num);
  pi_total.replace(datum::nan,1);
  return pi_total;
}

// [[Rcpp::export]]
arma::rowvec count_pi_double(arma::rowvec S, double i, NumericVector params, NumericVector constant) {
  
  double pi0 = params.containsElementNamed("gama0") ? params["gama0"] : constant["gama0"]; //containsElementNamed is good for our purposes
  double pi1 = params.containsElementNamed("gama1") ? params["gama1"] : constant["gama1"];
  double pi2 = params.containsElementNamed("gama2") ? params["gama2"] : constant["gama2"];
  
  
  arma::rowvec num =  exp(pi0 + pi1*i + pi2*exp(i - S));
  arma::rowvec pi_total =  num/(1+num);
  pi_total.replace(datum::nan,1);
  return pi_total;
}


// [[Rcpp::export]]
arma::mat count_p(double i, NumericVector params, NumericVector constant){
  arma::vec a = arma::linspace(1,i,i);
  arma::rowvec b = a.t();
  arma::mat S = arma::repmat(b,i,1);
  
  arma::mat thing = 1 - count_pi_mat(S,S.t(),params,constant);
  thing.row(i-1) =  arma::ones(i).t();
  arma::mat up = arma::ones(thing.n_rows,thing.n_cols);
  arma::mat uppertones =  triu(up,1);
  arma::mat thing2 = tril(thing,0) + uppertones;  
  arma::rowvec product = colProds(thing2);
  arma::rowvec pst = product % count_pi_double(b,i,params,constant);
  
  return(pst);
}

// [[Rcpp::export]]
arma::rowvec count_m(NumericVector i, NumericVector params, NumericVector constant){
  double m0 = params.containsElementNamed("beta0") ? params["beta0"] : constant["beta0"]; //containsElementNamed is good for our purposes
  double m1 = params.containsElementNamed("beta1") ? params["beta1"] : constant["beta1"];
  arma::rowvec m  = exp(m0 + m1*i);
  return(m);
}


// [[Rcpp::export]]
double count_log_like(NumericVector first_record_data,NumericVector params, NumericVector constant){
  
  NumericVector lambda (first_record_data.length());
  NumericVector summand2 (first_record_data.length());
  
  for(int i = 1; i < lambda.length() - 1; ++i){
    arma::rowvec S = arma::linspace(1,i,i).t();
    arma::rowvec Am = count_m(wrap(S),params,constant);
    arma::rowvec Ap = count_p(i,params,constant);
    lambda(i-1) = sum(Am % Ap);
    summand2(i-1) = first_record_data(i-1)*log(lambda(i-1)) - lambda(i-1);
  }
  double LL = -sum(summand2);
  return(LL);
}
