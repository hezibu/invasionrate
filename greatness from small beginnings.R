setwd("sac/")
Rcpp::compileAttributes(verbose=TRUE)

library(tools)
package_native_routine_registration_skeleton(dir = getwd())

RcppArmadillo::RcppArmadillo.package.skeleton("sacII")
Rcpp::compileAttributes(pkgdir = "sacII/",verbose=TRUE)

