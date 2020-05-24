function m = count_m(t,params)
%This function calculates the mean, mu, from Solow and Costello (2004)
m0 = params(1);
m1 = params(2);
m = exp(m0 + m1*t);
